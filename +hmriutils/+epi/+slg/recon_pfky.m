function [Irss yout] = recon_pfky(y, ny, method)
% function [Irss yout] = recon_pfky(y, ny, method)
%
% Partial Fourier reconstruction
%
% Inputs
%  y           [nx etl nz nc]     y(:,:,z,c) = 2D k-space for slice z, coil c
%  ny          [1]                full ky matrix size
%  method      string             'zerofill' (default) or... (TODO)
%
% Outputs
%  Irss        [nx ny nz]         root-sum-of-squares coil-combined image
%  yout        [nx ny nz nc]      k-space

if nargin < 3
    method = 'zerofill';
end

method = lower(method);

[nx etl nz nc] = size(y);

yout = zeros(nx, ny, nz, nc);

if strcmp(method, 'zerofill')
    yout(:,(end-etl+1):end,:,:) = y;
end

if strcmp(method, 'homodyne')
    for z = 1:nz
        for c = 1:nc
            yout(:,:,z,c) = hmriutils.epi.homodyne(y(:,:,z,c), ny, ...
                'trw_lp', 10, ...
                'trw_hp_lo', 6, ...
                'trw_hp_hi', 6);
        end
    end

end

Irss = zeros(nx,ny,nz);
for z = 1:nz
    [~, Irss(:,:,z)] = toppe.utils.ift3(squeeze(yout(:,:,z,:)), 'type', '2d');
end

