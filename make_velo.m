 pre_init
% if ~exist('change','var')
sinputs=0:elenum.S_input-1;
change1=rulecurr(sinputs+1)~=rulecurr(max(sinputs,1));
change0=rulecurr(sinputs+1)~=rulecurr(min(sinputs+2,18));
change2=rulecurr(sinputs+1)~=rulecurr(mod(sinputs+9,18)+1);
change=cat(1,change0,change1,change2)';
% change
% end
cells(x,y)=back;
% cells(x,y)=zeros(n,n);
% cells(2:2+pop-1,2:2+pop-1)=rand(pop,pop)<p0;
 % back=cells;
S_velo_old=zeros(n,n,2,lag+1);

k=k0;
rind=floor((k-1)/soupdensi)+1;    
getrule
middle=x(ceil(numel(x)/2));
S_defect=zeros(n,n);
S_defect(middle,middle)=1;
step=1;


subplot(2,2,2)
f3s=scatter3([],[],[],3,'x');

xlim([0 1])
ylim([0 n^2/10])
zlim([0 10])
%%
while step<10000
%%
step=step+1;
    cells=update(cells);

% S_popND=cat(3,ones(size(S_pop_old))*2,findneighbor(S_pop_old,neighbor.S_pop));
S_inputT=conv2(cells,FIR.S_input,'same');
S_input=S_inputT(x,y);
S_pop=cells(x,y);
S_popND=findneighbor(S_pop,neighbor.S_pop);
S_inputND=findneighbor(S_input,neighbor.S_pop);



make_change
%%
make_defect
% S_changeN=bsxfun(@rdivide,S_change,sum(S_change,3));
% S_changeN(isnan(S_changeN))=0;
% S_changeNND=findneighbor(S_changeN,neighbor.S_change);
% S_defectND=findneighbor(S_defect,neighbor.S_change);
% S_defect=sum(S_defectND.*S_changeNND,3);
area=sum(S_defect(:)~=0);
st=sum(S_defect(:));
pddx=(sum(S_defect,2)/st)';
pddy=sum(S_defect,1)/st;
meanx=sum(pddx.*(1:n));
meany=sum(pddy.*(1:n));
varx=sum(pddx.*((1:n)-meanx).^2);
vary=sum(pddy.*((1:n)-meany).^2);
tvar=varx+vary;
tstd=sqrt(tvar);
rec=[st,area,tstd];


subplot(2,2,1);
imagesc(S_defect')
subplot(2,2,2)
hold on
scatter3(st,area,tstd,3,'x');
hold off
%%
 
%rows=repmat(S_pop+1,1,1,9);
% cols=cat(3,S_input+18,S_inputND)+1;
% id2=sub2ind(siz,cols,rows);
S_flow=sum(S_changeinA,3);
% S_change(S_change==0)=nan;
S_velo=squeeze(sum(bsxfun(@times,S_change,shiftdim(neighbor.S_change,-2)),3,'omitnan'));
% S_velo=squeeze(mean(bsxfun(@times,S_change,shiftdim(neighbor.S_change,-2)),3,'omitnan'));
% S_velo=convn(S_velo,[1 1 1; 1 -9 1; 1 1 1]/9,'same');
% S_velo=convn(S_velo,ones(3,3),'same');
S_velo_old=cat(4,S_velo,S_velo_old(:,:,:,1:lag));
OS_velo=mean(S_velo_old,4);
S_speed=sum(abs(S_velo),3);




subplot(2,2,1);
imagesc(S_defect')
% imagesc((S_velo(:,:,1))')
% subplot(2,2,2);
% S_minspeed=S_speed==min(findneighbor(S_speed,neighbor.S_change),[],3);
% imagesc(S_minspeed')
% imagesc((S_velo(:,:,2))')
[x1,x2]=ndgrid(1:n,1:n);
subplot(2,2,3)
imagesc(S_pop')
% view(-90,90)
hold on 
quiver(x1,x2,OS_velo(:,:,1),OS_velo(:,:,2))

hold off
subplot(2,2,4)

imagesc(S_speed')

if mod(step,2)==0
drawnow
end
end
%%
S_change=bsxfun(@rdivide,S_change,sum(S_change,3));

% S_changeC=isnan(sum(S_change,3));
S_change(isnan(S_change))=1/9;
% S_change=cat(3,S_changeC,S_change);
S_changeND=findneighbor(S_change,neighbor.S_change);
% S_changeND_old=cat(4,S_changeND,S_changeND_old(:,:,:,1:lag));
% S_changeNDM=mean(S_changeND_old,4);
% S_defect=ones(n,n);
i=0;
while i<1;
S_defectND=findneighbor(S_defect,neighbor.S_change);
S_defect=sum(S_changeND.*S_defectND*mass,3)+S_defect*(1-mass);
i=i+1;
end
% S_defect=S_defect(2:end-1,2:end-1);
sum(S_defect(:));
D_defect=hist(S_defect(:),0:0.01:2);
