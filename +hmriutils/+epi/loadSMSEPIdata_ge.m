function dat = loadSMSEPIdata_ge(pfile, frame, etl, nz, mb)
%
% Load raw data acquired with smsepi.seq on GE scanners.
% Assumes that smsepi.seq is created with SMS-EPI/sequence/Pulseq/writeSMSEPI.m
%
% Output
%    dat    [nx etl npartitions ncoils]

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

% Get kspace locations and interpolate onto Cartesian grid
% Apply 1/2 sample shift to align echoes.
%system('tar xf smsepi.tar module3.mod');
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod('readout.mod');
kx = sysGE.raster*1e-6*sysGE.gamma*1e-4*cumsum(gx);  % cycles/cm
kx = kx - kx(end)/2;
kx = kx((hdr.npre+1):(hdr.npre+size(d,1))); % dwell time = 4us = gradient raster time, so number of samples = hdr.rfres here
del = 0.0;
kxo = interp1(1:nfid, kx, (1:nfid) - 0.5 - del, 'linear', 'extrap');
kxe = interp1(1:nfid, kx, (1:nfid) + 0.5 + del, 'linear', 'extrap');
kxe = fliplr(kxe);

fprintf('Interpolating... ');
fovXcm = fov(1)*100;
dco = hmriutils.epi.rampsamp2cart(d(:, 1:2:end, :, :), kxo, nx, fovXcm, 'nufft');   % odd echoes
dce = hmriutils.epi.rampsamp2cart(d(:, 2:2:end, :, :), kxe, nx, fovXcm, 'nufft');   % even echoes
fprintf(' done\n');
dat = zeros(nx, etl, npartitions, nCoils);
dat(:,1:2:end,:, :) = dco;
dat(:,2:2:end,:, :) = dce;
