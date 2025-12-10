function smask = getsamplingmask(ikz, nx, ny, mb,Ry)
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
%   Ry   [1]     in-plane undersample factor
%
% Output
%   smask    [nx ny mb]     3D k-space sampling mask (maskyz repeated at each kx point)

if nargin < 4
    mb = length(ikz);
    Ry=1;
elseif nargin < 5
    Ry = 1;
end

if length(ikz) == mb
    mb = length(ikz);
    tmp = [];
    for ii = 1:Ry:Ry*ny/mb
        tmp = [tmp ikz(:)'];
    end
    ikz = tmp;
else
    assert(~mod(ny, length(ikz)), 'ny must be a multiple of multiband factor (mb)');
end

smaskyz = false(Ry*ny,mb);
for y = 1:Ry:Ry*ny
    smaskyz(y, ikz(ceil(y/Ry))) = true;
end
smask = false(nx, Ry*ny, mb);
smask(:, smaskyz(:,1:mb)) = true;
