
S_popND=findneighbor(S_pop,neighbor.S_pop);

siz=size(change);
rows=repmat(S_input+1,1,1,9);
cols=cat(3,2*ones(n,n),S_popND)+1;
id=sub2ind(siz,rows,cols);

S_changeND=change(id);
S_change=findneighbor(S_changeND,-neighbor.S_change);
S_change=single(S_change);
% S_changeout=S_change(:,:,neighbor.project);
% S_changein=S_changeND;
% S_changeinA=S_changein-S_changeout;
