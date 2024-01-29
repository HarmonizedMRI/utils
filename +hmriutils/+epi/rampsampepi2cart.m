function dc = rampsampepi2cart(dr, kxo, kxe, nx, fovXcm, method)
% function dc = rampsampepi2cart(dr, kxo, kxe, nx, fovXcm, method)
%
% Interpolate ramp-sampled EPI to cartesian grid
%
%   dr       [nr, etl, ...]   Ramp-sampled raw data along 1st ("x") dimension, 
%                             EPI train along 2nd dimension.
%   kxo      [nr 1]           k-space sampling locations for odd echoes (cycles/cm)
%   kxe      [nr 1]           k-space sampling locations for even echoes (cycles/cm)
%   nx       [1]              Number of pixels (along x) in reconstructed image
%   fovXcm   [1]              FOV (cm) along x
%   method   string           'nufft' (default) or 'spline' (faster but less accurate)

import hmriutils.epi.*

if ~exist('method', 'var')
    method = 'nufft';
end

drSize = size(dr);
nr = drSize(1);
etl = drSize(2);

dr = reshape(dr, nr, etl, []);
dco = rampsamp2cart(dr(:,1:2:end,:), kxo, nx, fovXcm, method);   % odd echoes
dce = rampsamp2cart(dr(:,2:2:end,:), kxe, nx, fovXcm, method);   % even echoes
dc = zeros(nx, etl, size(dr,3));
dc(:,1:2:end,:) = dco;
dc(:,2:2:end,:) = dce;

dc = reshape(dc, [nx etl drSize(3:end)]);


