%%given
% S_input;
% SMI_vic;
%%
neighbor.S_flow=neighbor.S_pop;
% S_flow=ones(n,n);
S_inputND=findneighbor(S_input,neighbor.S_flow);
S_flowND=findneighbor(S_flow,neighbor.S_flow);
% nSMI_vic=SMI_vic/MI_vic;
% bonddict=logsig(SMI_vic);
bonddict=SMI_vic;
id=sub2ind(size(bonddict),S_inputND+1,repmat(S_input+1,1,1,length(neighbor.S_flow)));
bondND=bonddict(id);
bondNND=findneighbor(bondND,-neighbor.S_flow);
bondN=repmat(sum(bondNND,3),1,1,length(neighbor.S_flow));

lossexponent=-bondN(:,:,1);
HLE=std(lossexponent(:));

S_flowloss=2.^(-bondN(:,:,1));

bondN=bondNND./bondN.*(1-2.^(-bondN));
bond=findneighbor(bondN,-neighbor.S_flow);

% bond=b
% bond=bsxfun(@rdivide,bondND,sum(bondND,3));
% bond(:,:,1)=0;
bond=makezero(bond);

i=0;
% subplot(2,2,3)
% f3i=imagesc(S_flow);
% view(-90,90)

i=0;
while i<1
% set(f3i,'CData',S_flow)
% drawnow
% pause
S_flowND=findneighbor(S_flow,neighbor.S_flow);
% nS_flow=findneighbor(bsxfun(@times,bond,S_flow),-neighbor.S_flow);
nS_flow=bond.*S_flowND;
nS_flow=sum(nS_flow,3);
% dS_flow=sum(S_flowND.*bond,3);
dS_flow=nS_flow-S_flow;
% sum(abs(dS_flow(:))/900);
S_flow=nS_flow+1;
H_flow=log2(n^2)-sum(entropise(S_flow(:)));
% S_flow=S_flow+dS_flow;
sum(S_flow(:)); 
i=i+1;
end
D_flow=hist(S_flow(:),hedge);
LS_flow=S_flow(:)';
