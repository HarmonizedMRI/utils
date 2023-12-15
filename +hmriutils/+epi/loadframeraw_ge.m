function d = loadframeraw_ge(fn, etl, np, frame, acqOrder)
%
% Load SMS-EPI raw data for one temporal frame. 
% Assumes maxView = np*etl, unless acqOrder == true
%
% Inputs
%   fn               P-file name
%   etl              echo train length
%   np               number of partition/excitations (groups of SMS slices) per frame
%   frame            frame to load
%
% Input option
%   acqOrder         true or false (default)
%
% Ouput
%   d                [nFID etl np nc]

import hmriutils.io.*

if nargin < 5
    acqOrder = false;
end

if ~acqOrder
    echo = 1;
    sliceStart = frame + 1;
    sliceEnd = sliceStart;
    d = loadpfile(fn, echo, sliceStart, sliceEnd);    % [nfid ncoils 1 1 np*etl]
    d = flipdim(d, 1);                % tv6 flips data along FID direction
    d = squeeze(d);                   % [nfid nc np*etl]
    [nfid nc npxetl] = size(d);
    d = reshape(d, nfid, nc, etl, np);
    d = permute(d, [1 3 4 2]);   % [nfid etl np nc]
else
    d = toppe.utils.loadpfile(fn, [], [], [], 'acqOrder', true);    % [nfid ncoils nframes*np*etl]
    d = flipdim(d, 1);                % tv6 flips data along FID direction
    iStart = (frame-1)*np*etl + 1;
    iStop = iStart + np*etl-1;
    d = d(:, :, iStart:iStop);   % [nFID nc np*etl]
    [nfid nc N] = size(d);
    d = reshape(d, nfid, nc, etl, np);
    d = permute(d, [1 3 4 2]);   % [nfid etl np nc]
end

