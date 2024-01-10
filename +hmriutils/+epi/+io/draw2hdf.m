function draw2hdf(D, etl, np, ofn)
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

% sort data by frames
fprintf('Sorting data into %d frames...', nfr);
D = reshape(D, nfid, nc, etl, np, nfr);
D = permute(D, [1 3 4 2 5]);
fprintf(' done\n');

% Write to .h5 file. Delete if it already exists.
if isfile(ofn)
    delete(ofn)
end
fprintf('Writing .h5 file...');
h5create(ofn, '/kdata/real', [nfid etl np nc nfr], 'Datatype', class(D));
h5create(ofn, '/kdata/imag', [nfid etl np nc nfr], 'Datatype', class(D));
h5write(ofn, '/kdata/real', real(D));
h5write(ofn, '/kdata/imag', imag(D));
fprintf(' done\n');

return

% write one frame at a time. This fails in various ways -- Matlab hdf support seems flaky?
for ifr = 1:nfr
    d = hmriutils.epi.loadframeraw_ge(fn, etl, np, ifr, false);
    h5write(ofn, '/kdata/r', real(d), [1 1 1 1 ifr], [nfid etl np nc 1]);
%    h5write(ofn, '/kdata/i', imag(d), [1 1 1 1 ifr], [nfid etl np nc 1]);
end


