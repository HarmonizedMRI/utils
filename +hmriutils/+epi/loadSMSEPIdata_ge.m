function dat = loadSMSEPIdata_ge(pfile, frame, etl, nz, mb, readoutfile)
%
% Load ramp-sampled EPI data acquired with smsepi.seq on GE scanners,
% and interpolate onto cartesian grid.
% Assumes that smsepi.seq was created with SMS-EPI/sequence/Pulseq/writeSMSEPI.m
%
% Inputs
% ------
% pfile         File location of the P-file containing the SMS-EPI data
% frame         Temporal frame to read
% etl           echo train length
% nz            number of slices in reconstructed image volume
% mb            SMS/multiband factor
% readoutfile   Name of .mod file containing the EPI echo trapezoid waveform.
%
% Output
% ------
% data         k-space data for the given frame
%                Size: [nx etl npartitions ncoils], where npartitions = nz/mb, 
%                      i.e., the number of slices in the full reconstructed volume
%                      divided by the multiband factor
%
% Outputs
% -------
% data         k-space data for the given frame
%                Size: [nsamples, ncoils, npartitions], where npartitions is
%                      the number of slices in the full reconstructed volume
%                      divided by the multiband factor

% mask         k-space CAIPI sampling mask
%                Size: [nx, ny, mb] Note that this mask is the same for each set
%                      of simultaneously excited slices. Also note that
%                      nsamples = sum(mask(:))


% Load data for one frame
% nfid = number of ADC samples per echo (>nx due to ramp sampling)
% Note: nechoes = 1 by construction, see 'Pulseq on GE' manual
d = toppe.utils.loadpfile(fn, [], frame, frame);
d = flipdim(d, 1);      % [nfid ncoils 1 1 etl*npartitions]

d = permute(d, [1 2 5 3 4]);   % [nfid ncoils etl*npartitions] 

[nfid, ncoils, nviews] = size(d);

npartitions = nz/mb;
assert(nviews == etl*npartitions, ...
    'Number of views (%d) in data file does not equal etl*npartitions', nviews);

d = reshape(d, nfid, ncoils, etl, npartitions);

d = permute(d, [1 3 4 2]);   % [nfid etl npartitions ncoils]

% Get kspace locations (ramp sampled).
% 
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readoutfile);
gamma = 4.2576;      % kHz/Gauss
dt = 4e-3;   % GE gradient raster time, ms
kx = gamma*dt*cumsum(gx);   % cycles/cm
kx_ctr = mean(kx(end/2:(end/2)+1));
kx = kx - kx_ctr;    % assume that echo is sampled symmetrically about kx=0
kx = kx((hdr.npre+1):(hdr.npre+hdr.rfres));    % ADC on during this portion of the waveform
dur = numel(kx)*dt;
kxo = interp1([0 dt:dt:dur], [0 kx(:)'], dt/2:dt/2:dur);  % odd echoes
kxe = fliplr(kxo);

% Interpolate onto Cartesian grid
fprintf('Interpolating... ');
fovXcm = fov(1)*100;
dco = hmriutils.epi.rampsamp2cart(d(:, 1:2:end, :, :), kxo, nx, fovXcm, 'nufft');   % odd echoes
dce = hmriutils.epi.rampsamp2cart(d(:, 2:2:end, :, :), kxe, nx, fovXcm, 'nufft');   % even echoes
fprintf(' done\n');
dat = zeros(nx, etl, npartitions, nCoils);
dat(:,1:2:end,:, :) = dco;
dat(:,2:2:end,:, :) = dce;
