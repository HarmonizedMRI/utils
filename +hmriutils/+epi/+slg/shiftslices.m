function yout = shiftslices(yin, ph)
% Shift reconstructed slices away from/back to center of FOV
%
% Inputs:
%   yin     [nx ny mb nc]     Reconstructed images
%   ph      [ny mb]           Phase modulation along EPI train, see getsmsphasemodulation.m
%
% Output:
%   yout     [nx ny nz nc]

[nx ny mb nc] = size(yin);

for z = 1:mb
    for c = 1:nc
        for x = 1:nx
            yout(x,:,z,c) = yin(x,:,z,c).*exp(1i*ph(:,z).');
        end
    end
end
