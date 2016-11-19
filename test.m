
LS_inputo;

size(H_Second(LS_inputo(100,:),LS_inputo(1:100,:)))
dmat=squareform(pdist(LS_inputo,@(XI,XJ)H_Second(XI,XJ)));
