function Z_start = getsliceordering(np)
% Start slice order for the np partitions in SMS-EPI
% Interleaved, with last two shots/partitions swapped if np is even

Z_start = [1:2:np 2:2:np];
if ~mod(np,2)
    Z_start = [Z_start(1:end-2) Z_start(end) Z_start(end-1)]; 
end

