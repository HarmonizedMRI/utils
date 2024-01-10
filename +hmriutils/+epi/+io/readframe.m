function d = readframe(fn, frame)

sz = h5info(fn, '/kdata/real').Dataspace.Size;
dr = h5read(fn, '/kdata/real', [1 1 1 1 frame], [sz(1:4) 1]);
di = h5read(fn, '/kdata/imag', [1 1 1 1 frame], [sz(1:4) 1]);

d = double(dr) + 1i*double(di);
