

%% initialise cells
% figure('Visible','off')
    fig=gcf;
    fig.PaperPositionMode = 'auto';
    fig.Position=[0 0 800 800];
if ~exist('replayV','var')
    replayV=0;
end
global cells cells_old n x y rulecurr
n=40;
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

%%
middle=x(ceil(numel(x)/2));
s0.dfx=middle;
s0.dfy=middle;

particle=repmat(s0,partnum,1);
particle_old=repmat(particle,1,lag);
siz=size(particle_old);


fi=imagesc(cells);
tic
soup=1;
pdata_agg=zeros(lag,3,soupmax,'int8');

k=k0;


% sp0.pdata=[];
% sp=repmat(sp0,soupmax,1);
soupdensi=6;
k=k0-kstep;
kmax=size(ruletemp,1)*soupdensi;
kmax=1E100;
while k<kmax
    k=k+kstep;
% for k=k0:kstep:size(ruletemp,1)*soupdensi
      rind=floor((k-1)/soupdensi)+1;    
  
    getrule
    soup=1;
    
    while soup<=soupmax
if ~replayV
rsoup=single(rand(n,n)<p0);
end
cells(x,y)=rsoup;
particle=repmat(s0,partnum,1);
particle_old=repmat(particle,1,lag);

for step=1:cooldown
%     cells_old=cells;
    cells=update(cells);
end

% cells=zeros(n+2,n+2);
% cells(middle:middle+1,middle:middle+2)=[1 1 0; 1 1 1];

step=1;


while step<lag
    %%
    cells_old=cells;
    cellsdf=cells_old;
    cells=update(cells);
    step=step+1;    

    for pi=numel(particle):-1:1
        p=particle(pi);
        p=update_p(p);
        if p.ind==0
            particle(pi)=[];
            else
        particle(pi)=p;
        end
        cellsdf(p.dfx,p.dfy)=cellsdf(p.dfx,p.dfy)+5;
    end
    
    if ~isempty(particle)
    particle_old(1:size(particle,1),step)=particle;
    end
    set(fi,'CData',cellsdf)

    if true
    %     if toc>0.005
%     plot(mp)
%     if replayV
    if true
        drawnow
    pause(0.2)
    end
    tic
    end
    
end
    

xs=reshape([particle_old.dfx],siz);
ys=reshape([particle_old.dfy],siz);
ds0=cat(3,diff(xs,[],1),diff(ys,[],1));
ds=max(sum(abs(ds0),3),[],1);
mp0=cat(3,mean(xs,1),mean(ys,1))-middle;
mp=(sum(abs(mp0),3));

pdata_agg(:,:,soup)=[ds;mp;1:lag]';
soup
cond=ds(end)<7.5 && mp(end)>7.5;
if cond
        pause
end

% drawnow
soup=soup+1;

end



%% 
pdata=reshape(permute(pdata_agg,[2 1 3]),3,[]);
pdata=double(pdata);
ma=max(pdata(:));
mi=min(pdata(:));
po=double(ma+1);
count=zeros(ma+1,ma+1,ma+1);

l=size(pdata,2);
pdata=pdata+(rand(3,l)-0.5)/2;
xs=pdata(1,:);
ys=pdata(2,:);
zs=pdata(3,:);

subplot(2,2,3)
scatter3(xs,ys,zs,1.5,'x');
xlim([0 lag*2])
ylim([0 lag])
zlim([0 lag])
xlabel('diff')
ylabel('average pos')
zlabel('lag')
view(-90,0)

subplot(2,2,2)
scatter3(xs,ys,zs,1.5,'x');
xlim([0 lag*2])
ylim([0 lag])
zlim([0 lag])
xlabel('diff')
ylabel('average pos')
zlabel('lag')

view(0,0)
subplot(2,2,1)
scatter3(xs,ys,zs,1.5,'x');
xlim([0 lag*2])
ylim([0 lag])
zlim([0 lag])
xlabel('diff')
ylabel('average pos')
zlabel('lag')

view(-90,90)
subplot(2,2,4)
scatter3(xs,ys,zs,1.5,'x');
xlim([0 lag*2])
ylim([0 lag])
zlim([0 lag])
xlabel('diff')
ylabel('average pos')
zlabel('lag')

view(-30,30)


% suptitle(['defect particle mobility' ' rulenumber='  num2str(rind) ' ' rulename{rind}])
makeprint
end
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