
% LS_defectL=log2(S_defect(:)');
% %     S_change=rule_change(S_popND);
%     change=rule(S_inputP+1)~=rule(S_input+1);
%     S_defect=change;
%     S_defectT(x,y)=S_defect;
%     S_defectT(1,:)=S_defectT(x(end),:);
%     S_defectT(end,:)=S_defectT(x(1),:);
%     S_defectT(:,1)=S_defectT(:,y(end));
%     S_defectT(:,end)=S_defectT(:,y(1));
% sinputs=0:elenum.S_input-1;
% change1=rulecurr(sinputs+1)~=rulecurr(max(sinputs,1));
% change0=rulecurr(sinputs+1)~=rulecurr(min(sinputs+2,18));
% change2=rulecurr(sinputs+1)~=rulecurr(mod(sinputs+9,18)+1);
% change=cat(1,change0,change1,change2)';

%%
% S_popND=cat(3,ones(size(S_pop_old))*2,findneighbor(S_pop_old,neighbor.S_pop));
% S_popND=findneighbor(S_pop_old,neighbor.S_pop);


% S_popND=S_popND(2:end-1,2:end-1,:);
% siz=size(change);
% rows=repmat(padarray(S_input,[1 1],0,'both'),1,1,9)+1;
% rows=repmat(S_input+1,1,1,9);
% cols=cat(3,2*ones(n,n),S_popND)+1;
% cols=S_popND+1;
% id=sub2ind(siz,rows,cols);
% S_changeND=change(id);
% S_change=findneighbor(S_changeND,-neighbor.S_change);
% S_change=cat(3,9-sum(S_change,3),S_change);
% S_change=bsxfun(@rdivide,S_change,sum(S_change,3));
S_changeN=bsxfun(@rdivide,S_change,sum(S_change,3));
S_changeN(isnan(S_changeN))=0;
if ~exist('mass','var')
    mass=1
end
% S_changeC=isnan(sum(S_change,3));
% S_change(isnan(S_change))=1/9;
% S_change=cat(3,S_changeC,S_change);
S_changeNND=findneighbor(S_changeN,neighbor.S_change);
% S_changeND_old=cat(4,S_changeND,S_changeND_old(:,:,:,1:lag));
% S_changeNDM=mean(S_changeND_old,4);
% S_defect=ones(n,n);
i=0;
while i<1;
S_defectND=findneighbor(S_defect,neighbor.S_change);
S_defect=sum(S_changeNND.*S_defectND*mass,3)+S_defect*(1-mass);
i=i+1;
end
% S_defect=S_defect(2:end-1,2:end-1);
% sum(S_defect(:));
% D_defect=hist(S_defect(:),0:0.01:2);
