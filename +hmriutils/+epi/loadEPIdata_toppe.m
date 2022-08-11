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

nz = size(dat, 2);
nCoils = size(dat,3);
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

% crop turns
nCrop = floor(epiInfo.nBlipMax/2);
kx = kx((nCrop+1):(end-nCrop),:);
dat = dat((nCrop+1):(end-nCrop), :, :);

% interpolate onto Cartesian grid
kxo = kx(:,1);
kxe = kx(:,2);
fprintf('Interpolating... ');
dco = rampsamp2cart(dat(:, 1:2:end, :), kxo, nx, fov);   % odd echoes
dce = rampsamp2cart(dat(:, 2:2:end, :), kxe, nx, fov);   % even echoes
fprintf(' done\n');
dat = zeros(nx, ny, nz*nCoils);
dat(:,1:2:end,:) = dco;
dat(:,2:2:end,:) = dce;

% Reshape
dat = reshape(dat, [nx ny nz nCoils]);

return

