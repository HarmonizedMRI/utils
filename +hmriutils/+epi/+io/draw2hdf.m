function draw2hdf(D, etl, np, ofn, varargin)
%
% Write SMS-EPI fMRI raw data matrix to .h5 file,
% after reshaping by (temporal) frame.
%
% Inputs
%   D      [nfid nc N]      SMS/3D EPI raw data matrix (any type), for N sequential RF shots/FIDs.
%                           nfid = number of data points per shot/FID.
%                           nc = number of receive coils
%                           Number of temporal frames is N/(etl*np).
%   etl                     echo train length
%   np                      number of partition/excitations (groups of SMS slices) per frame
%   ofn                     .h5 output file name
%
% Input options (keyboard arguments)
%   maxFramesPerFile        Split data set into multiple .h5 files. 
%                           This seems to be needed for large D, perhaps due to a bug in h5write()
% 
% Output
%   ofn   .h5 file containing the data matrix, reshaped to size [nfid etl np nc nframes]
%
% Example: Load a GE P-file, write data to .h5 file, and read a frame from it: 
%   >> D = toppe.utils.loadpfile('P12345.7', [], [], [], 'acqOrder', true, 'returnAsDouble', false);
%   >> hmriutils.epi.io.draw2hdf(D, 72, 10, 'mydata.h5');
%   >> d = hmriutils.epi.io.readframe('mydata.h5', 47);   % size(d) = [nfid etl np nc]

[nfid, nc, N] = size(D);
nfr = N/etl/np;    % number of temporal frames
assert(~mod(nfr,1), 'size(D,3) must be multiple of etl*np');

arg.maxFramesPerFile = nfr;
arg = toppe.utils.vararg_pair(arg, varargin);

fnStem = hmriutils.epi.io.getfilenamestem(ofn);

% Exit if file already exists
if isfile([fnStem '.h5'])
    fprintf('File already exists -- exiting\n');
    return
end

% sort data by frames
fprintf('Sorting data into %d frames...', nfr);
D = reshape(D, nfid, nc, etl, np, nfr);
D = permute(D, [1 3 4 2 5]);
fprintf(' done\n');

% Write header/entry .h5 file
nFiles = ceil(nfr/arg.maxFramesPerFile);
ofn = [fnStem '.h5'];

fprintf('Writing %s...', ofn);
h5create(ofn, '/maxFramesPerFile', [1], 'Datatype', class(arg.maxFramesPerFile));
h5write(ofn, '/maxFramesPerFile', arg.maxFramesPerFile);

h5create(ofn, '/nFiles', [1], 'Datatype', class(nFiles));
h5write(ofn, '/nFiles', nFiles);

nFramesLastFile = nfr - floor(nfr/arg.maxFramesPerFile)*arg.maxFramesPerFile;
h5create(ofn, '/nFramesLastFile', [1], 'Datatype', class(nFramesLastFile));
h5write(ofn, '/nFramesLastFile', nFramesLastFile);

h5create(ofn, '/dataSize', [ndims(D)], 'Datatype', class(size(D)));
h5write(ofn, '/dataSize', size(D));
fprintf(' done\n');

% write data file(s)
for ii = 1:nFiles
    ofn = [fnStem '_' num2str(ii) '.h5'];

    % frames to write to file
    FR = (ii-1)*arg.maxFramesPerFile+1 : ii*arg.maxFramesPerFile;
    if ii == nFiles
        FR = FR(1:nFramesLastFile);
    end

    % write to file
    fprintf('Writing %s...', ofn);
    h5create(ofn, '/kdata/real', [nfid etl np nc length(FR)], 'Datatype', class(D));
    h5create(ofn, '/kdata/imag', [nfid etl np nc length(FR)], 'Datatype', class(D));
    h5write(ofn, '/kdata/real', real(D(:,:,:,:,FR)));
    h5write(ofn, '/kdata/imag', imag(D(:,:,:,:,FR)));
    fprintf(' done\n');
end

return

% write one frame at a time. This fails in various ways -- Matlab hdf support seems flaky?
for ifr = 1:nfr
    d = hmriutils.epi.loadframeraw_ge(fn, etl, np, ifr, false);
    h5write(ofn, '/kdata/r', real(d), [1 1 1 1 ifr], [nfid etl np nc 1]);
%    h5write(ofn, '/kdata/i', imag(d), [1 1 1 1 ifr], [nfid etl np nc 1]);
end


