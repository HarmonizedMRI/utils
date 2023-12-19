function [Irss, w] = recon(ysms, ycal, Z, nz, smask, K, w)
% function [Irss, w] = recon(ysms, ycal, Z, nz, smask, K, [w])
%
% Reconstruct SMS data with split slice GRAPPA
%
% Inputs:
%   ysms    [nx ny nc]        SMS EPI data (one RF shot)
%   ycal    [nx ny mb nc]     Single-slice, unshifted 'ACS' k-space data for slice GRAPPA calibration.
%                             Non-zero values define the calibration region.
%   Z       [mb]              SMS slice indices in full image volume
%   nz      [1]               Number of slices in full image volume
%   smask   [nx ny mb]        Sampled (3D) k-space locations along echo train,
%                             which defines the z blips. See getsamplingmask.m
%   K       [2]               Kernel size (e.g., [5 5] or [7 7])
%
% Optional input:
%   w       {mb nc} cell array    w{z,c} = slice GRAPPA weights for slice z, coil c.
%                                 If not provided, calculate and return to caller.
%
% Outputs:
%   Irss    [nx ny mb]        Root-sum-of-squares coil-combined reconstructed image
%   w                         slice GRAPPA weights

[nx ny mb nc] = size(ycal);

% phase-correct sms data for slice groups not centered at iso-center
ysms = hmriutils.epi.slg.smsphasecorrect(ysms, Z(1), nz, smask);

% apply FOV shift to calibration data to match SMS acquisition
ph = hmriutils.epi.slg.getsmsphasemodulation(nz, Z, smask);
ycal = hmriutils.epi.slg.shiftslices(ycal, ph);

% calculate slice GRAPPA weights (if not provided by caller)
if ~exist('w', 'var')
    tic
    [~, pinvA] = hmriutils.epi.slg.slgcal(ycal, 1, 1, K, 'lam', 1.^2);  % get pinv(A)
    for z = 1:mb
        fprintf('.');
        for c = 1:nc
            w{z,c} = hmriutils.epi.slg.slgcal(ycal, z, c, K, 'pinvA', pinvA);
        end
    end
    fprintf('\nCalibration time: %.2fs\n', toc);
end

% reconstruct
tic
y = zeros(size(ycal));
for z = 1:mb
    for c = 1:nc
        % reconstruct slice z, coil c
        for c2 = 1:nc
            y(:,:,z,c) = y(:,:,z,c) + conv2(ysms(:,:,c2), w{z,c}(:,:,c2), 'same');
        end
    end
end
fprintf('Recon time: %.2fs\n', toc);

% shift reconstructed slices back to center of FOV
y = hmriutils.epi.slg.shiftslices(y, -ph);

% Return root-sum-of-squares coil-combined image
Irss = zeros(nx, ny, mb);
for z = 1:mb
    [~, Irss(:,:,z)] = toppe.utils.ift3(squeeze(y(:,:,z,:)), 'type', '2d');
end

