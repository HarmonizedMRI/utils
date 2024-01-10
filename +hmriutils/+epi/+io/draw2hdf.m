function draw2hdf(draw, etl, np, ofn)
%
% Write SMS-EPI fMRI raw data matrix to .h5 file.
%
% Inputs
%   draw    [nfid nc N]     Raw data matrix (single)
%   etl                     echo train length
%   np                      number of partition/excitations (groups of SMS slices) per frame
%   ofn                     .h5 output file name
%
% Example loading GE P-file: 
%   >> draw = toppe.utils.loadpfile('P12345.7', [], [], [], 'acqOrder', true, 'returnAsDouble', false);
%   >> hmriutils.epi.io.draw2hdf(draw, 72, 10, 'mydata.h5');

[nfid, nc, N] = size(draw);

nfr = N/etl/np;    % number of temporal frames

assert(~mod(nfr,1), 'size(draw,3) must be multiple of etl*np');

% delete .h5 file if already existing
if isfile(ofn)
    delete(ofn)
end

% sort data into frames
fprintf('Sorting data into %d frames...', nfr);
D = single(zeros(nfid, etl, np, nc, nfr));   % full data matrix
for ifr = 1:nfr
    % d = data for one temporal frame
    iStart = (ifr-1)*np*etl + 1;
    iStop = iStart + np*etl-1;
    d = draw(:, :, iStart:iStop);   % [nFID nc np*etl]
    [nfid nc N] = size(d);
    d = reshape(d, nfid, nc, etl, np);
    d = permute(d, [1 3 4 2]);   % [nfid etl np nc]. Data
    D(:,:,:,:,ifr) = d;
end

% write to .h5 file
if isfile(ofn)
    delete(ofn)
end
fprintf('\nWriting .h5 file...');
h5create(ofn, '/kdata/real', [nfid etl np nc nfr])
h5create(ofn, '/kdata/imag', [nfid etl np nc nfr])
h5write(ofn, '/kdata/real', real(D));
h5write(ofn, '/kdata/imag', imag(D));
fprintf('done\n');

return

% write one frame at a time. This fails in various ways -- Matlab hdf support seems flaky?
for ifr = 1:nfr
    d = hmriutils.epi.loadframeraw_ge(fn, etl, np, ifr, false);
    h5write(ofn, '/kdata/r', real(d), [1 1 1 1 ifr], [nfid etl np nc 1]);
%    h5write(ofn, '/kdata/i', imag(d), [1 1 1 1 ifr], [nfid etl np nc 1]);
end


