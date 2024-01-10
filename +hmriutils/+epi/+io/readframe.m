function d = readframe(fn, frame)
%
% Load raw data for one frame from .h5 file created with draw2hdf.m.
% Return as double.
%
% Inputs
%   fn      string   .h5 file name, see draw2hdf
%   frame   int       frame number to load
%
% Output
%   d       [nfid etl np nc]    k-space (raw) data for one frame

sz = h5info(fn, '/kdata/real').Dataspace.Size;
dr = h5read(fn, '/kdata/real', [1 1 1 1 frame], [sz(1:4) 1]);
di = h5read(fn, '/kdata/imag', [1 1 1 1 frame], [sz(1:4) 1]);

d = double(dr) + 1i*double(di);
