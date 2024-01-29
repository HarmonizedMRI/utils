function d = readframe(fn, frame)
%
% Load raw data for one frame from .h5 file created with draw2hdf.m.
% Return as double.
%
% Inputs
%   fn      string   .h5 file name (with/without '.h5' extension), see also draw2hdf
%   frame   int       frame number to load
%
% Output
%   d       [nfid etl np nc]    k-space (raw) data for one temporal frame

maxFramesPerFile = h5read(fn, '/maxFramesPerFile');
nFiles = h5read(fn, '/nFiles');
dataSize = h5read(fn, '/dataSize');
dataSize = dataSize(:)';   % must be row vector

d = zeros(dataSize(1:end-1));  % size of one temporal frame

% get name of file containing requested frame
fnStem = hmriutils.epi.io.getfilenamestem(fn);
iFile = ceil(frame/maxFramesPerFile); 
ifn = [fnStem '_' num2str(iFile) '.h5'];

% read one frame
ifr = frame - (iFile-1)*maxFramesPerFile;
fprintf('Reading frame %d from %s\n', ifr, ifn); 
dr = h5read(ifn, '/kdata/real', [1 1 1 1 ifr], [dataSize(1:4) 1]);
di = h5read(ifn, '/kdata/imag', [1 1 1 1 ifr], [dataSize(1:4) 1]);
d = double(dr) + 1i*double(di);

return


%% old code
for ii = 1:nFiles
    ofn = [fnStem '_' num2str(ii) '.h5'];
    FR = (ii-1)*arg.maxFramesPerFile+1 : ii*arg.maxFramesPerFile;
    if ii == nFiles
        FR = FR(1:nFramesLastFile);
    end
    fprintf('Reading %s...', ofn);

end

sz = h5info(ifn, '/kdata/real').Dataspace.Size;
dr = h5read(ifn, '/kdata/real', [1 1 1 1 frame], [sz(1:4) 1]);
di = h5read(ifn, '/kdata/imag', [1 1 1 1 frame], [sz(1:4) 1]);
