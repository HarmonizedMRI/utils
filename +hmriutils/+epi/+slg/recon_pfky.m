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
%  Irss        [nx ny nz nc]      root-sum-of-squares coil-combined image
%  yout        [nx ny nz nc]      k-space

if nargin < 3
    method = 'zerofill';
end

method = lower(method);

[nx etl nz nc] = size(y);

if strcmp(method, 'zerofill')
    yout = zeros(nx, ny, nz, nc);
    yout(:,(end-etl+1):end,:,:) = y;
end

Irss = zeros(nx,ny,nz);
for z = 1:nz
    [~, Irss(:,:,z)] = toppe.utils.ift3(squeeze(yout(:,:,z,:)), 'type', '2d');
end

