

%% initialise cells
% figure('Visible','off')
    fig=gcf;
    fig.PaperPositionMode = 'auto';
    fig.Position=[0 0 800 800];
if ~exist('replayV','var')
    replayV=0;
end
global cells cells_old n x y rulecurr async
n=50;
partnum=100;
x=2:n+1;
y=2:n+1;
cells=zeros(n+2,n+2);

% soup=single(rand(n,n)<p0);
% cells(x,y)=soup;

if ~exist('cooldown','var')
cooldown=10;
end

%%
middle=x(ceil(numel(x)/2));

subplot(2,2,1)
fi=imagesc(cells);

subplot(2,2,2)
fh=histogram2([],[],0:9,0:9);
zlim([0 200]);

subplot(2,2,3)
f3i=imagesc(cells);
caxis([0 9])
subplot(2,2,4)
f4i=imagesc(cells);
caxis([0 9])

tic
soup=1;
pdata_agg=zeros(soupmax,3);
% pdata_agg=zeros(lag,3,soupmax,'int8');

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
preallocate
    soup=1;
    
while soup<=soupmax
    if ~replayV
        rsoup=single(rand(n,n)<p0);
    end
    cells=zeros(n+2,n+2);
    cells(middle:middle+pop,middle:middle+pop)=rsoup(1:pop+1,1:pop+1);
%     cells(x,y)=rsoup;
    S_defect=zeros(n,n);
    S_defect(middle,middle)=1;


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
        step=step+1;    
        S_pop=cells(x,y);
        S_inputT=conv2(cells,FIR.S_input,'same');
        S_input=S_inputT(x,y);
        cells=rulecurr(S_inputT+1);
        cells=torus(cells);

        make_change
        make_defect
        S_changeout=sum(S_change(:,:,neighbor.project),3);
        S_changein=sum(S_changeND,3);
        
%         sum(S_defect(:))
        set(f3i,'CData',S_changein)
        set(f4i,'CData',S_changeout)

%         set(fh,'Data',[S_changein(:),S_changeout(:)])
        if true
        %     if toc>0.005
    %     plot(mp)
%         if replayV
        if mod(step,2)==0
            drawnow
        pause(0.2)
        end
        tic
        end

    end
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
    massd=abs(meanx-middle)+abs(meany-middle);
    rec=[st,massd,tstd];

%     pdata(soup,:)=rec;
    pdata_agg(soup,:)=rec;
    soup
    cond=ds(end)<7.5 && mp(end)>7.5;
    if cond
        pause
    end

    % drawnow
    soup=soup+1;

    end



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