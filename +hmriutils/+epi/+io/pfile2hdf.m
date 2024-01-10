function pfile2hdf(fn, etl, np, nfr, ofn)
%
% Load pfile acquired with HarmonizedMRI SMS-EPI sequence,
% and write data to .h5 file.
%
% Inputs
%   pfile            P-file name
%   etl              echo train length
%   np               number of partition/excitations (groups of SMS slices) per frame
%   nfr              total number of frames in acquisition
%   ofn              .h5 output file name
%
% Example: 

% get data size
d = hmriutils.epi.loadframeraw_ge(fn, etl, np, 1, false);

[nfid, etl_, np_, nc] = size(d);

assert(etl_ == etl, 'etl provided does not match P-file');
assert(np_ == np, 'num partitions provided does not match P-file');

% delete file if already existing
if isfile(ofn)
    delete(ofn)
end

% initialize data arrays
h5create(ofn, '/kdata/r', [nfid etl np nc nfr])
%h5create(ofn, '/kdata/i', [nfid, etl, np, nc, nfr])

% write one frame at a time
for ifr = 1:1
    d = hmriutils.epi.loadframeraw_ge(fn, etl, np, ifr, false);
    h5write(ofn, '/kdata/r', real(d), [1 1 1 1 ifr], [nfid etl np nc 1]);
%    h5write(ofn, '/kdata/i', imag(d), [1 1 1 1 ifr], [nfid etl np nc 1]);
end

return

% create output file
h5create(ofn, '/kdata_r', size(kdata))
h5create(ofn, '/kdata_i', size(kdata))

% write data to file
disp('Writing kspace data to file...')
h5write(arg.outfile, '/kdata_r', real(kdata));
h5write(arg.outfile, '/kdata_i', imag(kdata));
disp('Done.')

