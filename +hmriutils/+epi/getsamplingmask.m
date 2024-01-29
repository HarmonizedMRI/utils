function smask = getsamplingmask(ikz, nx, ny, mb)
% Define SMS EPI sampling mask
%
% Input
%   ikz    [mb] or [ny]   kz sampling indices along echo train.
%                         If length(ikz) == mb, ikz is repeated ny/mb times.
%   nx     [1]            matrix size
%   ny     [1]            matrix size
%
% Optional input
%   mb   [1]     SMS multiband factor. If not provided, assume mb = length(ikz)
%
% Output
%   smask    [nx ny mb]     3D k-space sampling mask (maskyz repeated at each kx point)

if nargin < 4
    mb = length(ikz);
end

if length(ikz) == mb
    mb = length(ikz);
    tmp = [];
    for ii = 1:ny/mb
        tmp = [tmp ikz(:)'];
    end
    ikz = tmp;
else
    assert(~mod(ny, length(ikz)), 'ny must be a multiple of multiband factor (mb)');
end

smaskyz = false(ny,mb);
for y = 1:ny
    smaskyz(y, ikz(y)) = true;
end
smask = false(nx, ny, mb);
smask(:, smaskyz(:,1:mb)) = true;
