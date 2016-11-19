if ~exist('replayV','var')
    replayV=1;
end
replayV=1-replayV
eval(vmethod)