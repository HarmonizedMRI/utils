function d = loadframe_ge(fn, etl, np, frame)
%
% Load SMS-EPI raw data for one temporal frame. Assumes maxView = np*etl.
%
% Inputs
%   fn               P-file name
%   etl              echo train length
%   np               number of partition/excitations (groups of SMS slices) per frame
%   frame            frame to load
%
% Ouput
%   d                [nFID etl np nc]

import hmriutils.io.*

echo = 1;
sliceStart = frame + 1;
sliceEnd = sliceStart;
d = loadpfile(fn, echo, sliceStart, sliceEnd);    % [nfid ncoils 1 1 np*etl]
d = flipdim(d, 1);                % tv6 flips data along FID direction
d = squeeze(d);                   % [nfid nc np*etl]
[nfid nc npxetl] = size(d);
d = reshape(d, nfid, nc, etl, np);
d = permute(d, [1 3 4 2]);   % [nfid etl np nc]

