function dout = smsphasecorrect(d, z_start, nz, smask)
% function dout = smsphasecorrect(d, z_start, nz, smask)
%
% Phase-correct SMS EPI data. This is necessary when the
% center (mb/2+1) slice in the SMS group is offset 
% from the scanner iso-center.
% 
% Inputs
%   d          [nx ny nc]      One CAIPI EPI shot
%   nz         [1]             Number of slices in full image volume
%   z_start    [1]             Slice index of first slice in SMS slice group (partition)
%   smask      [nx ny mb]      k-space sampling mask (view EPI shot as undersampled 3D k-space data)
%
% Output
%   dout      [nx ny nc]

[nx ny mb] = size(smask);
nc = size(d,3);

% Prepare the data for partition (SMS slice group) p.
% With no CAIPI sampling scheme, all the k-space data for a partition lies
% in a single kz plane. However, with CAIPI sampling, each ky line
% actually comes from a different kz plane. To account for this,
% first we replicate the k-space data, ending up with mb identical
% kz planes. Then we use the sampling mask to keep only the k-space
% samples that were actually acquired in each kz plane.
d3d = zeros(nx, ny, mb, nc);
for c = 1:nc
    d3d(:,:,:,c) = repmat(d(:,:,c), 1, 1, mb);
end
for c = 1:nc
    d3d(:,:,:,c) = d3d(:,:,:,c).*smask;
end

np = nz/mb;                         % number of partitions (slice groups)
z_iso = nz/2 + 1;                   % center of full volume
z_p = z_start + np * floor(mb/2);   % center of excited slice group (partition)
z_offset = z_p - z_iso;             % offset from isocenter

deltak = 1/nz;     % kz-space spacing, in units of 1/pixel

for ii = 1:mb
    kz = (-mb/2+ii-1)*deltak;   % kz level, units of 1/pixel
    d3d(:,:,ii,:) = d3d(:,:,ii,:) * exp(1i*2*pi*kz*z_offset);
end

dout = squeeze(sum(d3d,3));
