function dat = loadEPIdata_toppe(pfile, frame, epiInfo, sys)
% function dat = loadEPIdata_toppe(pfile, frame, epiInfo, sys)
%
% Load ramp-sampled EPI data and return data on Cartesian grid.
% Assumes sequence created with caipiepifmri.m
%
% Inputs:
%   pfile      GE raw data file
%   frame      time frame to load
%   epiInfo    struct containing info needed to load and reconstruct
%              See caipiepifmri.m
%   sys        TOPPE system struct
%
% Output:
%   dat   [nx ny nz nCoils]

N = epiInfo.N;
FOV = epiInfo.FOV;
gpre = epiInfo.gpre;   % x prewinder (G/cm)
gx = epiInfo.gx;       % full echo train (G/cm)
gx1 = epiInfo.gx1;     % one echo

[dat, rdb_hdr] = toppe.utils.loadpfile(pfile);
dat = flipdim(dat, 1);  
echo = 1;
dat = squeeze(dat(:,:,:,echo,frame));   % [nFID nCoils nz]
dat = permute(dat, [1 3 2]);            % [nFID nz nCoils]

nz = size(dat, 2)
nCoils = size(dat,3)
nFID = size(dat,1);

% get kspace locations
% apply 1/2 sample shift to align echoes
nx = N(1); ny = N(2); fov = FOV(1);
kx = sys.raster*sys.gamma*(cumsum(gx) - sum(gpre.x));  % cycles/cm
kx = interp1(1:nFID, kx, (1:nFID) - 0.5);

% sort into 2d arrays
nr = length(gx1);
dat = dat(1:(nr*ny), :);
kx = kx(1:(nr*ny));
dat = reshape(dat, nr, ny, []);
kx = reshape(kx, nr, ny);

% fix first k-space location (NaN due to interp1)
kx(1,:) = kx(2,:);

% interpolate onto Cartesian grid
kxo = kx(:,1);
kxe = kx(:,2);
dco = rampsamp2cart(dat(:, 1:2:end, :), kxo, nx, fov);   % odd echoes
dce = rampsamp2cart(dat(:, 2:2:end, :), kxe, nx, fov);   % even echoes
dat = zeros(nx, ny, nz*nCoils);
dat(:,1:2:end,:) = dco;
dat(:,2:2:end,:) = dce;

return

%% Do odd/even echo correction

% reconstruct echoes and interpolate to Cartesian grid
% first get system nufft object A and density function (for speed)
addpath ~/github/HarmonizedMRI/SMS-EPI/recon/ % reconecho.m
[~,Ao,dcfo] = reconecho([], nx, [], [], kx(:,1), fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kx(:,2), fov); % even echoes
x = zeros(nx,ny,nCoils); % data after 1D NUFFT along kx
dCart = zeros(nx, ny, nCoils);  % Cartesian data
for ic = 1:1:nCoils
	for iy = 1:2:ny
		x(:,iy,ic) = reconecho(dat(:,iy,ic), nx, Ao, dcfo);
	end
	for iy = 2:2:ny
		x(:,iy,ic) = reconecho(dat(:,iy,ic), nx, Ae, dcfe);
	end
    for iy = 1:ny
	    dCart(:,iy,ic) = fftshift(fft(fftshift(x(:,iy,ic),1), [], 1),1);
    end
end

%th = zeros(nx,1);  % odd/even echo phase difference
% th = th + abs(tmp).^2.*exp(1i*angle(xtmp(:,2:2:,ic)./xtmp(

%% Do odd/even echo correction

% Estimate off-resonance phase from even echoes (center line)
off = zeros(1,ny); % radians
for ic = 1:nCoils
    x1 = x(end/2,2:2:end,ic);
    tmp = angle(circshift(x1,1)./x1);
%    off = abs(x1(1:length(off))).^2*exp(1i*tmp);
    size(off)
%    off = mean(off(2:(end-2)));
%    off = angle(exp(1i*x(end/2,2:2:end,ic));
end

return

pho = angle(x(end/2,1:2:end));
dph = angle(exp(1i*phe)./exp(1i*pho)); % / (2*pi*n*sys.raster);
offresPhase = mean(dph(2:end)); % exclude first echo
x(:,1:2:end) = x(:,1:2:end).*exp(1i*offresPhase);
%dCart(:,1:2:end,:) = dCart(:,1:2:end,:).*exp(1i*phOffset);

% get odd/even echo linear phase offset
phoe = angle(x(:,2:2:end)./x(:,1:2:end));
tmp = angle(x2./x1.*exp(1i*angle(x1./x3)/2));
th = th + abs(x1).^2.*exp(1i*tmp);
xsos = xsos + abs(x1).^2;

th = zeros(nx,1);
xsos = zeros(nx,1);  % sum-of-squares coil combined image (for mask)
	
for ic = 1:1:nCoils
    x1 = x(:,1:2:end);
    x2 = x(:,2:2:end);

	tmp = angle(x2./x1.*exp(1i*angle(x1./x3)/2));
	th = th + abs(x1).^2.*exp(1i*tmp);
	xsos = xsos + abs(x1).^2;
end

th = angle(th);
xsos = sqrt(xsos);

mask = xsos > 0.1*max(xsos(:));

% fit phase difference to 2d plane
H = [ones(sum(mask(:)),1) X(mask)];  % spatial basis matrix (2d linear)
ph(isl,:) = H\th(mask);  

thhat = embed(H*ph(isl,:)', mask);
figure; plot([th thhat th-thhat]); legend


