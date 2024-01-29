function fnStem = getfilenamestem(fn)
% get part of file name without the '.h5' extension if present

if strcmp(fn(end-2:end), '.h5')
    fnStem = fn(1:end-3);
else
    fnStem = ofn;
end
