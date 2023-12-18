function ph = getsmsphasemodulation(nz, Z, smask)
%
% Get phase modulation due to gz blips along the EPI train.
% This is used to shift reconstructed slices back to center of FOV
% (or more generally, to undo effect of z blips).
%
% Inputs
%   nz      [1]             Number of slices in full image volume
%   Z       [mb]            SMS slice indices in full image volume
%   smask   [nx ny mb]      Sampled (3D) k-space locations along echo train,
%                           which defines the z blips. See getsamplingmask.m
% Output
%   ph      [ny mb]         Phase modulation

[nx ny mb] = size(smask);

assert(length(Z) == mb, 'Number of slices (Z) and size(smask,3) must match');

% 'Excite' one slice at a time and apply z blips as defined in 'smask',
% and get resulting phase modulation 
d_ex_1 = zeros(nx, ny, mb);
for z = 1:mb
    d_ex_1(:,:,z) = 1;  % 'excite' a slice
    [~, tmp] = hmriutils.epi.slg.simulatesmsepishot(d_ex_1, nz, Z, smask);
    ph(:,z) = tmp(:,z);
    d_ex_1(:,:,z) = 0;  % undo to get ready for next slice
end

