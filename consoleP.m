%%

tmax=500;
tol1=1/10;
amax=1;
%%
Rd=10;
R=ceil(Rd+tmax*(amax+tol1+1));
t0=1;
%%
tota=0.1;
aseq=-amax:tota:amax;
[ax ay]=ndgrid(aseq,aseq);
% as=[ax(:) ay(:)];
as=cat(3,ax,ay);
[a1,a2]=size(ax);

dsetx0=R-Rd+1:R+Rd;
dsety0=dsetx0;
dsetmin=min(dsetx0);
dsetmax=max(dsetx0);

[abx ]=ndgrid([dsetmin dsetmax]);
% ab=cat(3,abx,aby);
% ab=shiftdim(reshape(ab,[],2)',-2);
ab=shiftdim(abx,-3);
[uabx ]=ndgrid([-1 1]);
% uab=cat(3,uabx,uaby);
% uab=shiftdim(reshape(uab,[],2)',-2);
uab=shiftdim(uabx,-3);
% ep=reshape(200*tol1);
%%
S_pop=zeros(2*R,2*R);
S_input=zeros(2*R,2*R);
S_defect=zeros(2*R,2*R);
S_defect(dsetx0,dsety0)=1;
maxLEA=zeros(a1,a2);
% S_defect=2*R
l1=2*R;
l2=2*R;
xset=1:l1;
yset=1:l2;
S_pop=rand(size(S_pop))<0.5;
figure(4)
subplot(2,2,1);
f1i=imagesc(S_pop);
subplot(2,2,2)
f2i=imagesc(maxLEa);
subplot(2,2,3)
f3i=imagesc(S_pop);
r=importrule({'B012/S1256'});
rulecurr=[r{1},r{2}];
t=t0;
preallocate

LLEA=[];

while t<=tmax
xsetn=xset(2:end-1);
ysetn=yset(2:end-1);
xshift=xsetn(1)-1;
yshift=ysetn(1)-1;
% S_input=conv2(single(S_pop(xset,yset)),FIR.S_input,'valid');
% S_pop=rulecurr(S_input(xsetn,ysetn)+1);
S_pop_old=S_pop;
S_input_old=S_input;
S_input=conv2(single(S_pop_old),FIR.S_input,'valid');
S_pop=rulecurr(S_input+1);

% S_defect=conv2(single(S_defect),FIR.S_input,'valid');


S_popND=cat(3,ones(size(S_pop_old))*2,findneighbor(S_pop_old,neighbor.S_pop));
% S_popND=S_popND(2:end-1,2:end-1,:);
siz=size(change);
rows=repmat(padarray(S_input,[1 1],0,'both'),1,1,9)+1;
cols=S_popND+1;
id=sub2ind(siz,rows,cols);
S_change=change(id);
S_defectND=cat(3,S_defect,findneighbor(S_defect,neighbor.S_pop));
S_defect=sum(S_change.*S_defectND,3);
S_defect=S_defect(2:end-1,2:end-1);


%     % S_change=bsxfun(@(sinput,state) change(sinput+1,state+1),S_input,S_popND);
% S_popND=S_popND(2:end-1,2:end-1);
% change0(xsetn,ysetn)=rulecurr(S_input(xsetn,ysetn)+1)~=rulecurr(max(S_input(xsetn,ysetn),1));
% change1(xsetn,ysetn)=rulecurr(S_input(xsetn,ysetn)+1)~=rulecurr(min(S_input(xsetn,ysetn)+2,18));
% change01=cat(4,change0,change1);
% %   cchange=rulecurr(S_input+1)~=rulecurr(mod(S_input+9,18)+1);
%     
%     cchange=rulecurr(1:elenum.S_input)~=rulecurr(mod((0:elenum.S_input-1)+9,18)+1);
%     S_popNDV(xsetn,ysetn,:,:)=bsxfun(@eq,S_popND(xsetn,ysetn,:),reshape(0:1,1,1,1,2));
%     S_change(xsetn,ysetn,:)=sum(bsxfun(@times,S_popNDV(xsetn,ysetn,:,:),change01(xsetn,ysetn,:,:)),4);
%     S_defectND(xset,yset,:)=findneighbor(S_defect(xset,yset),neighbor.S_pop);
%     S_defect(xsetn,ysetn,:)=sum(cat(3,cchange(S_input(xsetn,ysetn)+1),S_change(xsetn,ysetn,:)).*cat(3,S_defect(xsetn,ysetn),S_defectND(xsetn,ysetn,:)),3);
%     tic
%     S_pop(xsetn,ysetn)=rulecurr(S_input(xsetn,ysetn)+1);
%     t1=toc;
%     tic
%     a=S_input(xsetn,ysetn);
%     b=rulecurr(a+1);
%     toc-t1
% S_pop=rulecurr(S_input+1);
acone=bsxfun(@plus,bsxfun(@plus,as*t,ab),uab*floor(tol1*t));
squeeze(acone(11,11,:,:))
c1=ceil(acone(:,:,1,1));
c2=floor(acone(:,:,1,2));
c3=ceil(acone(:,:,2,1));
c4=floor(acone(:,:,2,2));
LEA=arrayfun(@(a,b,c,d) makezero(log2(sum(reshape(S_defect((a:b)-xshift,(c:d)-yshift),1,[]))))/t,c1,c2,c3,c4);
maxLEA=max(maxLEA,LEA);
t=t+1;
xset=xsetn;
yset=ysetn;
disp(sprintf('t=%d',t))
LLEA=[LLEA;LEA(:)'];
set(f1i,'CData',S_defect>0);
set(f2i,'CData',LEA);
% set(f3i,'CData',sum(S_change,3));
% 
drawnow
pause(0.02)
% universe_act=universe(xset,yset);
end
%%
figure(5)
maxLEA=reshape(max(LLEA,[],1),21,21);
subplot(2,2,2)
% surface(LEA)
imagesc(LEA)
axc=gca;
axc.XTickLabels=num2str(aseq(axc.XTick)');
axc.YTickLabels=num2str(aseq(axc.YTick)');

subplot(2,2,[3 4])
imagesc(LLEA);
