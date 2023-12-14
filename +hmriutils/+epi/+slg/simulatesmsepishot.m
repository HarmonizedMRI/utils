function [ysms, ph] = simulatesmsepishot(d_ex, nz, IZ, smask)
%
% Simulate EPI sampling (one shot) of SMS group 'd_ex'.
%
% Inputs
%   d_ex    [nx ny mb nc]    Excited slice group (k-space) after excitation
%   nz      [1]              Number of slices in full image volume
%   IZ      [mb]             Slice indices in full image volume
%   smask   [nx ny mb]       sampled k-space locations (view readout as sampling 3D k-space),
%                            see getsamplingmask.m
%
% Output
%   ysms    [nx ny nc]       Simulated data for one EPI train
%   ph      [ny mb]          Phase difference between input and output.
%                            This is used to shift reconstructed slices back to center of FOV
%                            (or more generally, to undo effect of z blips), 
%                            see also getsmsphasemodulation.m

[nx ny mb nc] = size(d_ex);

d_ex = permute(d_ex, [1 2 4 3]);

deltak = 1/nz;    % in units of 1/pixel (ok since voxel size divides out in phase modulation step)
iz_iso = nz/2+1;  % center of slice group 

% loop over kz encoding levels within slice group
d_ex_blipped = zeros(size(d_ex));
for ikz_g = -mb/2:mb/2-1  

    % apply phase modulation created by z blip
    kz_g = ikz_g*deltak;      % kz-space encoding level in unit of 1/pixel
    dtmp = zeros(size(d_ex));
    for iz_g = 1:mb             
        kz = kz_g*(IZ(iz_g)-iz_iso); 
        dtmp(:,:,:,iz_g) = d_ex(:,:,:,iz_g) * exp(-1i*2*pi*kz); 
    end

    % collapse slices
    d_ex_blipped(:,:,:,ikz_g+mb/2+1) = sum(dtmp, 4);
end

d_ex_blipped = permute(d_ex_blipped, [1 2 4 3]);   % [nx ny mb nc]

dsms3d = 0*d_ex_blipped;
for c = 1:nc
    dsms3d(:,:,:,c) = d_ex_blipped(:,:,:,c) .* smask;
end
% dsms3d now represents the SMS-EPI data for one shot, as a 3D matrix.
% size(dsms3d) = [nx ny mb nc], with zeros where not sampled

% collapse slices
ysms = squeeze(sum(dsms3d, 3));   % [nx ny nc]

% return phase difference w.r.t. input, after phase correction
% (feels like there's an easier way to get ph but this works)
ph = zeros(ny, mb);
ysmspc = hmriutils.epi.slg.smsphasecorrect(ysms, IZ(1), nz, smask);
c = 1;
for z = 1:mb
    ph(:,z) = angle(ysmspc(nx/2,:,c)./d_ex(nx/2,:,c,z));
end

return
