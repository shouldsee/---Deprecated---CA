function[p]=update_p(p)        
global cells_old cells  x y n
cellsdftd=cells_old;                              
cellsdftd(p.dfx,p.dfy)=1-cells_old(p.dfx,p.dfy);

%          cellsdftd=update(cellsdftd,rulecurr);
dfds=xor(cells,update(cellsdftd));
ind_dfds=find(dfds(x,y));
if isempty(ind_dfds)
    p.ind=0;
    return
    % do nothing
else
    ind=datasample(ind_dfds,1);
    [dfx,dfy]=ind2sub([n n],ind);
    dfx=dfx+1;
    dfy=dfy+1;
%     p.ds=p.ds+abs(dfx-p.dfx)+abs(dfy-p.dfy);
    p.dfx=dfx;
    p.dfy=dfy;
end

%     p.pt=pt;
%     p.ddfx=p.dfx-p.dfx0;
%     p.ddfy=p.dfy-p.dfy0;
end