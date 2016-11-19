
global cells cells_old n x y rulecurr async
rind=1;
k=rind;
randrule=1;
randrule=0;
% name='nerd';
loadrule
getrule;
n=200;
x=2:n+1;
y=2:n+1;
cells=zeros(n+2,n+2);
async=0;
siz=size(cells);
[x1,x2]=ndgrid(x,y);
xyid=sub2ind(siz,x1,x2);
p0=0.5;

% n=5;
%%
rec=[];
rind=1;
loadrule;
n=4;
while rind<size(rulename,1);
% ind_all=(1:size(lndict,3));
% inputs=squeeze(convn(ndict,FIR.S_input,'valid'));
getrule
rsvct=2.^((1:(n-2)^2)-1);
imax=2^(n^2);
lndict=zeros(3,imax);
i0=1;
%%
T=zeros(2*rsvct(end),2*rsvct(end));
fir=reshape(rsvct,n-2,n-2);
parfor i=i0:imax
    OtS=reshape(str2num(dec2base(i-1,2,n^2)')',n,n);
%     OtSc=OtS(2:end-1,2:end-1);
    num=conv2(OtS,fir,'valid');
    OtS=rulecurr(1+conv2(OtS,FIR.S_input,'valid'));
%     OtSc=conv2(OtS,FIR.S_input,'same')
%     num=sum(OtSc(:).*rsvct(:));
    
%     conv2(OtSc,FIR.S_input,'valid');
%     Snum=sum(OtS(:).*rsvct(:));
  Snum=conv2(OtS,fir,'valid');

    dnum=hist(num(:),0:(rsvct(end)*2-1));
    dSnum=hist(Snum(:),0:(rsvct(end)*2-1));
    ddS=dnum'*dSnum;
    ddS=ddS/sum(ddS(:));
    %     
    T=T+ddS;
%     lndict(:,i)=[num,Snum,num+2*rsvct(end)*Snum];
    fprintf('%d of %d \n',i,imax)
end
%%
% d=hist(lndict',(0:(rsvct(end)*2)^2-1));
% d=bsxfun(@rdivide,d,imax*ones(1,3));
% HD=-log2(d).*d;
% T=reshape(d(:,3),2*rsvct(end),2*rsvct(end));
% H=sum(HD,1,'omitnan');
% H(4)=H(1)+H(2)-H(3);
eg=eig(T);
rec(rind,:)=eg;
fprintf('maxeig=%s, rind=%d, %s \n',num2str(abs(eg(1)),3),rind,rulename{rind});
rind=rind+1;
end
%%
fig=figure(3);
loadrule;
% hold on 
% scatter(rec(:,2),rec(:,3));
plot(sum(abs(ak.(name)),2));
dcm_obj=datacursormode(fig);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,rulename})
hold off

ak.(name)=rec;

%%
fn=fieldnames(ak);
figure(3)
cla
rules=[];
recs=[];
for i=1:size(fn,1);
   rec=ak.(fn{i});
    recs=[recs;rec];
   hold on
   scatter(rec(:,2),rec(:,3));
   name=fn{i};
    loadrule;
    rules=[rules;rulename];
end
figure(4)
clf
fig=gcf;
scatter(recs(:,2),recs(:,3));
dcm_obj=datacursormode(fig);

set(dcm_obj,'UpdateFcn',{@myupdatefcn,rules})
   
