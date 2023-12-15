function IZ_start = getsliceordering(np)
% Start slice order for the np partitions in SMS-EPI
% Interleaved, with last two shots/partitions swapped

IZ_start = [1:2:np 2:2:np];
IZ_start = [IZ_start(1:end-2) IZ_start(end) IZ_start(end-1)]; 

