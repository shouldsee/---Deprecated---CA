

%% initialise cells
% figure('Visible','off')
    fig=gcf;
    fig.PaperPositionMode = 'auto';
    fig.Position=[0 0 800 800];
if ~exist('replayV','var')
    replayV=0;
end
global cells cells_old n x y rulecurr
if ~exist('n','var')
n=3;
end
partnum=100;
x=2:n+1;
y=2:n+1;
cells=zeros(n+2,n+2);
% soup=single(rand(n,n)<p0);
% cells(x,y)=soup;

clear s0
s0.dfx=1;
s0.dfy=1;

s0.ind=1;
if ~exist('cooldown','var')
cooldown=10;
end


% middle=x(ceil(numel(x)/2));
s0.dfx=middle;
s0.dfy=middle;



subplot(2,2,1)
% fi=imagesc(cells);
fi=imagesc([]);
subplot(2,2,2)
fh=histogram2([],[],0:9,0:9);
zlim([0 200]);
subplot(2,2,3)
f3i=imagesc(cells);
caxis([0 9])
subplot(2,2,4)
f4i=imagesc(cells);
% caxis([0 9])
tic



soupnum=1;

rulemax=size(ruletemp,1);
rind=0;
rulemax=1E100;
siz=size(cells);
[x1,x2]=ndgrid(x,y);
xyid=sub2ind(siz,x1,x2);
binsoup=zeros(n^4,n^2);
dd=zeros((n^(2*2)),1);

ind2soup=@(a) str2num(dec2base(a,2,n^2)')';
soup2ind=@(a) base2dec(num2str(reshape(a,[],1))',2);

for soup=1:2^(n^2)
    
bsoup=ind2soup(soup-1);
binsoup(soup,:)=bsoup;
end
sdmat=squareform(pdist(binsoup,'hamm'));
[s1,s2]=ndgrid(1:2^(n^2),1:2^(n^2));
% s1=s1(:);
% s2=s2(:);

rind=rind0;
%%
while rind<rulemax
% for k=k0:kstep:size(ruletemp,1)*soupdensi
    rind=rind+1;
    rind=rind;
%     rind=floor((k-1)/soupdensi)+1;    
  
    getrule
    preallocate
    soupnum=1;
    %%
    traj_all=zeros(2^(n^2),2^(n^2));
for soup=1:2^(n^2)
    traj=zeros(2^(n^2),1);
    %%initiate soup
    % binsoup(soup,:)=bsoup;
    soupind=soup;
    cells(xyid)=binsoup(soup,:);
    while traj(soupind)==0;
    
    traj(soupind)=1;
    cells_old=cells;
    cells=update(cells);
    soupind=soup2ind(cells(xyid))+1;
    end
    traj_all(soup,:)=traj;
%     df=cells_old==cells;
%     dpri=mean(mean(df(xyid)));
%     dd(soup)=min(dpri,1-dpri);    

end

%%
traj_sum=sum(traj_all,1);
subplot(2,2,1)
plot(traj_sum)

subplot(2,2,2)
histogram(repmat(traj_sum,2,1));
ylim([1 2^(n^2)*2])
set(gca,'YScale','log')
subplot(2,2,3)
histogram(sum(traj_all,2));
k=rind;
suptitle(sprintf('Rule %s, sum(all traj)=%d',rulename{rind},sum(traj_sum)));
drawnow

% makeprint
pause
%%

end


%% hamm dist visualiser
% set(fi,'CData',dmat)
% caxis(fi.Parent,[0 0.5])
% xlim(fi.Parent,[0 n^2])
% ylim(fi.Parent,[0 n^2])
% 
% subplot(2,2,2)
% xs=sdmat(:);
% zs=xs+(rand(size(xs))-0.5)*0.05;
% xs=dd(s1(:))+(rand(size(xs))-0.5)/20;
% ys=dd(s2(:))+(rand(size(xs))-0.5)/20;
% 
% 
% scatter3(xs,ys,zs,5,'x');
% scatter(zs,abs(xs-ys),5,'x');
% xlim([-0.1 0.6])
% ylim([-0.1 0.6])
% zlim([-0.2 1.2])
% subplot(2,2,3)
% histogram2(dd(s1(:)),dd(s2(:)));
% zlim([1 (2^(n^2))^2]);
% 
% set(gca,'ZScale','log')
% subplot(2,2,4)
% M=(abs(dd(s1)-dd(s2)))+0*sdmat;
% [V D]=eigs(M,2,'SA');
% [V,ind]=sort(V);
% 
% set(f4i,'CData',M(ind,ind));
% 
% xlim(f4i.Parent,[0 2^(n^2)])
% ylim(f4i.Parent,[0 2^(n^2)])
% caxis(f4i.Parent,[0 0.5])


    %% 
%     pdata=reshape(permute(pdata_agg,[2 1 3]),3,[]);
%     pdata=double(pdata);
%     ma=max(pdata(:));
%     mi=min(pdata(:));
%     po=double(ma+1);
%     count=zeros(ma+1,ma+1,ma+1);
    pdata=pdata_agg';
    l=size(pdata,2);
%     pdata=pdata+(rand(3,l)-0.5)/5;
    xs=pdata(1,:);
    ys=pdata(2,:);
    zs=pdata(3,:);
    xlab='survived defect';
    ylab='displacement of mass centre';
    zlab='total std';
    xlimit=[-0.1 1.1];
    ylimit=[-1 lag*2];
    zlimit=[-1 lag];
    subplot(2,2,3)
    scatter3(xs,ys,zs,1.5,'x');
    xlim(xlimit)
    ylim(ylimit)
    zlim(zlimit)
    xlabel(xlab)
    ylabel(ylab)
    zlabel(zlab)
    view(-90,0)

    subplot(2,2,2)
    scatter3(xs,ys,zs,1.5,'x');
    xlim(xlimit)
    ylim(ylimit)
    zlim(zlimit)
    xlabel(xlab)
    ylabel(ylab)
    zlabel(zlab)

    view(0,0)
    subplot(2,2,1)
    scatter3(xs,ys,zs,1.5,'x');
    xlim(xlimit)
    ylim(ylimit)
    zlim(zlimit)
    xlabel(xlab)
    ylabel(ylab)
    zlabel(zlab)

    view(-90,90)
    subplot(2,2,4)
    scatter3(xs,ys,zs,1.5,'x');
    xlim(xlimit)
    ylim(ylimit)
    zlim(zlimit)
    xlabel(xlab)
    ylabel(ylab)
    zlabel(zlab)

    view(-30,30)


    % suptitle(['defect particle mobility' ' rulenumber='  num2str(rind) ' ' rulename{rind}])
    makeprint

%%
% for i=1:size(pdata,2)
% 
%     
% end
    %%
% hpdata=bsxfun(@times,double(pdata),permute(po.^(0:2),[2 1]));
% shpdata=sum(hpdata,1);
% % hhpdata=[;hist(shpdata,0:max(shpdata))];
% 
% %%
% 
% cord=mod(min(shpdata):max(shpdata),po^2);
% hpdata=bsxfun(@eq,pdata,shiftdim(mi:ma,-1));
% mhpdata=bsxfun(@times,hpdata(1,:,:),permute(hpdata(2,:,:),[1 2 4 3]));
% mhpdata=bsxfun(@times,mhpdata,permute(hpdata(3,:,:),[1 2 5 4 3]));
% mhpdata=squeeze(mean(hpdata,2));