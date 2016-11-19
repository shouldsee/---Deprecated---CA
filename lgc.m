if ~exist('stepnum0','var')
    stepnum0=1;
end

rec=[];

bd_act=0;

tic
lag=20;
global gifcount figs
while rind<=size(rulename,1);zlim([0 0.4])
cells=zeros(n+2,n+2);

n3=ceil(n/3);
cells(n3+1:2*n3,n3+1:2*n3)=single(rand(n3,n3)<0.55);
cold=repmat(cells,1,1,lag);
rdold=repmat(cells(xyid),1,1,lag2);
gh=fspecial('gaussian');
edf=[-1 -1 -1; -1 8 -1; -1 -1 -1]/8;
% edf=fspecial('laplacian',0)
sigmas=linspace(0.1,2,2);
egs=[];
% loadrule
% ruletemp=importrule(rulename);
getrule;
stepmax=1000;
figure(1)
subplot(4,1,1)
% la=plot(1:length(sigmas));
heg=linspace(0,0.2,100);
hs=hist(cells(:),heg);
% figs.h=plot(hs);
figs.h=histogram([0 25],100);
set(gca,'YScale','log')
ylim([1 1E4])
subplot(4,1,2:4)
figs.fi=imagesc(cells);
figure(2)
% clf
% figs.fi1=imshow(cells,[-0,0.2]);
figs.fi1=imagesc(cells);
colormap gray
% figs.fi1=surface(sgrd(xyid));
% view(40,70)
% zlim([0 0.4])
figure(3)
clf
figs.fi2=imagesc(cells,[-0 1]*0.2);
% figs.fi2=imagesc(cells);
gifcount=1;
title(figs.fi2.Parent,rulename{rind});

figure(5)
figs.h3=histogram2([0 0.2],[0 30],100);
view(60,60)

% ylim([0 1]/20)

% rulename{rind}
sgrd_old=repmat(cells,1,1,lag);
sgrd=cells;
md=floor(lag/2);
opticalFlow =opticalFlowHS;
figure(4)
colormap gray
for stepnum=stepnum0:stepmax;
    S_input=conv2(cells,FIR.S_input,'same');
    cells=rulecurr(S_input+1);
    cells=torus(cells);
    
    if stepnum >cooldown && mod(stepnum,6)==0;
    gcls=imgaussfilt(cells, sigmas(end));
    sgrd_old=cat(3,sgrd,sgrd_old(:,:,1:lag-1));
    sgrd=stdfilt(gcls);       
%     hs=hist(sgrd(:),heg);
%     hso(stepnum,:)=hs;    
%     
    sgrdm=mean(sgrd_old(:,:,1:6),3);   
    flow=estimateFlow(opticalFlow,sgrdm);
    agflow=conv2(flow.Magnitude,FIRM,'same')./(rangefilt(flow.Vx)+rangefilt(flow.Vy));
    dsgrdm=imhmax(sgrdm,0.05);
    dsgrdmmax=max(findneighbor(dsgrdm,neighbor.S_change),[],3);
    fagflow=( dsgrdmmax<0.0003).*agflow;   
    fagflow(isnan(fagflow))=0;gcls=imgaussfilt(cells, sigmas(end));
    sgrd_old=cat(3,sgrd,sgrd_old(:,:,1:lag-1));
    sgrd=stdfilt(gcls); 
    bd_act=mean(sgrd(border));
    
%     hs=hist(sgrd(:),heg);
%     hso(stepnum,:)=hs;    
%     
    sgrdm=mean(sgrd_old(:,:,1:6),3);   
%     flow=estimateFlow(opticalFlow,sgrdm);
%     agflow=conv2(flow.Magnitude,FIRM,'same')./(rangefilt(flow.Vx)+rangefilt(flow.Vy));
    dsgrd=imhmax(sgrd,0.05);
   
    randmat=rand(n+2,n+2);
%     if sum(cells(oxid))>numel(oxid)/6;
    if    mod(stepnum,360)==0
    cells(oxid)=((dsgrd(oxid)<1E-5)+((dsgrd(oxid)>1E-5).*randmat(oxid)<0.0)) & cells(oxid);
    end
    if mod(stepnum,360)==0;
    cells(ixid)=randmat(ixid)<0.5;
    end
    
%     cells=torus(cells);
cells=absorb(cells);
%     dsgrdm=imhmax(sgrdm,0.05);

%     dsgrdmmax=max(findneighbor(dsgrdm,neighbor.S_change),[],3);
%     fagflow=( dsgrdmmax<0.0003).*agflow;   
%     fagflow(isnan(fagflow))=0;
%     imagesc(sgrdm)
%     hold on
%         plot(flow,'DecimationFactor',[5 5],'ScaleFactor',5*2500)
%     hold off
%     agflow=conv2(flow.Magnitude,FIRM,'same')./(sqrt(stdfilt(flow.Vy)^2+stdfilt(flow.Vx)^2));
%     set(figs.fi1,'CData',(sgrdm-imhmax(sgrdm,0.05,4)))    
    %     agflow=flow.Magnitude./(rangefilt(flow.Vx)+rangefilt(flow.Vy));
%      agflow=flow.Magnitude./(sqrt(stdfilt(flow.Vy)^2+stdfilt(flow.Vx)^2));
%     bw0=sgrdm/max(sgrdm(:));
%     bw=sgrdm>0.03;
%     bw2=imfill(bw,'holes');
%     bw2_perim=bwperim(bw2);
%     bw_perim=bwperim(bw);
%     overlay1 = imoverlay(bw0, bw2_perim, [.3 1 .3]);
%     bw=bw0/max(bw0(:))>0.2;
%     unit=ones(n+2,n+2);
%     area=splitapply(@sum,unit(:),L(:)+1);
%     vxs=splitapply(@mean,flow.Vx(:),L(:)+1);
%     vys=splitapply(@mean,flow.Vy(:),L(:)+1);
%     avs=abs(vxs)+abs(vys);
%     imagesc(overlay1)
%     set(figs.fi1,'ZData',sgrdm(xyid));
%       set(figs.fi1,'CData',label2rgb(L));
%     imagesc(imextendedmax(sgrdm,0.1))
%     fagflow_old=fagflow;
%     dsgrdm=imhmax(sgrdm,0.05);
%     dsgrdmmax=max(findneighbor(dsgrdm,neighbor.S_change),[],3);
% %     fagflow=(sgrdm-imhmax(sgrdm,0.075,4)>0 & dsgrdm<1E-3).*agflow;
%     fagflow=( dsgrdmmax<0.0003).*agflow;

%     ofagflow=imopen(fagflow,strel('disk',3));
%     fagflow_peri=bwperim(fagflow);
%     L=watershed(ofagflow);
    set(figs.fi1,'CData',(sgrd));   
    set(figs.fi2,'CData',dsgrd>0);
%     set(figs.fi2,'CData',label2rgb(watershed(dsgrd)));
%     set(figs.fi2,'CData',fagflow);
%     set(figs.fi2,'CData',fagflow-fagflow_old,ones(10,10)));
%     set(figs.h,'Data',max(fagflow-5,0));

    drawnow
%     writegif('fi1',sprintf('./gallery/edges/%s/%s',name, rulename{rind}));
%      writegif('fi2',rulename{rind});

     pause(0.01)
   end
    
%     pause(0.01)
    %     cold=cat(3,cells,cold(:,:,1:lag-1));
    %     rdold=cat(3,rd,rdold(:,:,1:lag2-1));
%     set(fi,'CData',diff(rdold(:,:,1:2),[],3)/2^pw);
end
    rs_cold=convn(cold,rsvct,'valid');
    rrs_cold=squeeze(reshape(rs_cold,1,(n+2)^2,[]));
    d=sum(hist(rrs_cold',0:2^pw-1)'~=0,2);
    rd=reshape(d,n+2,n+2);
    nrd=rd/2^pw;
%     nrd=nrd-conv2(nrd,ones(3,3)/9,'same');
%     nrd=torus(nrd);
%     nrd=bsxfun(@minus,nrd,min(nrd(:)));
%         set(figs.fi2,'CData',S_input);
    for si=1:length(sigmas)
        grd=imgaussfilt(nrd, sigmas(si));
        egrd=conv2(grd,edf,'same');
        egrd=torus(egrd);
        gcls=imgaussfilt(cells, sigmas(si));
        sgrd=stdfilt(gcls);
       
%         egs(si)=mean(abs(egrd(:)));
        egs(si)=mean(abs(sgrd(:)));
%         S_egs=(sgrd>0.03)&(egrd<1E-5);
        S_egs=sgrd;
        hs=hist(S_egs(:),heg);

%          egs=mean(S_egs(:));

        S_egs=torus(S_egs);
        set(figs.fi,'CData',(sgrd));
        set(figs.fi1,'CData',S_egs);
        set(figs.fi2,'CData',abs(egrd));
%         set(figs.h,'YData',hs);
        
        drawnow
        
        pause(0.1)
    end
    S_input=S_input(xyid);
    ds=hist(S_input(:),0:17);
    H=sum(entropise(ds));
%     [~,si]=min(abs(egs));
%         grd=imgaussfilt(nrd, sigmas(si));
%         egrd=conv2(grd,edf,'valid');
%         egs(si)=mean(abs(egrd(:)));
    degs=abs([0 diff(egs)]);
    pwer=sigmas*degs';
%     set(la,'XData',sigmas);
%     title(sprintf('power=%.4f',pwer))
    mnrd=mean(nrd(:));
%     rec(rind,:)=[egs,-100,H,mnrd,pwer];
%     rec(rind,:,:)=shiftdim(hso,-1);
%     rec(rind,:)=egs;
    fprintf('egs=%s, rind=%d, %s \n',num2str(egs),rind,rulename{rind});
%     pause
    
%     set(fi,'CData',egrd);
%     drawnow
% pause

rind=rind+1;
%%
if mod(rind,10)==0;
    toc;
    tic;

end

end

%%
ak.(name)=rec;
if strcmp(name,'rand');
    save('rand','rulename');
end
rs.(name)=rulename;



%%
% names={'nerd','weirdo','rand'};
names={'nerd'};
% names={'nerd','rand'};
rules=[];
recs=[];
group=[];

for i=1:numel(names);
    name=names{i};
    rec=ak.(name);
    recs=[recs;rec];
    rulename=rs.(name);
    rulename=rulename(1:size(rec,1));
    rules=[rules;rulename];
    group=[group;i*ones(size(rec,1),1)];
		
end
figure(2)
[coeff,score,latent,tsquared,explained,mu1]=pca(squeeze(reshape(recs,1,[],length(heg))));
rscore=reshape(score(:,1:2),[],cooldown,2);
ind=find(mean(mean(abs(diff(rscore(:,end-5:end,:),[],2)),2),3)>1E1);
rscore=rscore(ind,:,:);
rules=rules(ind);
group=group(ind);

t=repmat(1:cooldown,size(rscore,1),1,2);
group=repmat(group,1,cooldown,2);
rules=repmat(rules,1,cooldown,2);
% (recs-mu1)/coeff;

cmap=colormap('lines');
xs=rscore(:,:,1);
ys=rscore(:,:,2);
zs=t(:,:,1);
gs=group(:,:,1);
scatter3(xs(:),ys(:),zs(:),15,cmap(gs(:),:));
% xlim([0.1,4E4])
% scatter3(score(:,1),score(:,2),score(:,3),15,cmap(group,:));
% gscatter(score(:,1),score(:,2),group);
fig=gcf;
dcm_obj=datacursormode(fig);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,rules(:)});

% save('agg','ak');
% save('agg','rs');
% %%
% load('agg','ak');









