
n=30;
x=2:n+1;
y=2:n+1;
eso=[];
egs=[];
v1s=[];
v2s=[];
timp_D_inputo=[];
timp_H_popo=[];
H_ego=[];
inds=[];
E_dmat_V_inputos=[];
lists={'v1','v2'};
scalars={'D_input','ANFlux','NFlux','RFlux','H_input',+...
     'H_pop','D_Tinput','VSD_Tinput','V_input','v1','v2','v3'};
% cmap=colormap('colorcube');
if ~exist('dirname')
    dirname='gallery/temp/';
end
mkdir(dirname)
fprintf('saveing to %s \n',dirname);
trigger=1;
if ~exist('k0')
    k0=1;
end

for k=k0:length(S_poplist)
    for i=1:length(scalars)
        eval([scalars{i} 'o=[];'])
       
    end
    lag=2;
    axs=[];
    inds=[];
    S_popo=S_poplist{k};
    if sum(S_popo(:))==0
        for i=1:length(lists)
           list=lists{i};
           eval(sprintf('%ss=[%ss;zeros(1,size(%ss,2))];',list,list,list))
        end
    
        continue
    end

S_popo=reshape(S_popo,[],30,30);
cells=zeros(n+2,n+2);
D_transio=[];
D_inputo=[];
for stepnumber=1:size(S_popo,1)
    fprintf('\n clip k=%d,step %d ',k,stepnumber)
    S_pop=squeeze(S_popo(stepnumber,:,:));
    cells(x,y)=S_pop;
    cells(1,:)=cells(x(end),:);
    cells(end,:)=cells(x(1),:);
    cells(:,1)=cells(:,y(end));
    cells(:,end)=cells(:,y(1));
    
    D_pop=hist(S_pop(:),0:1);
    H_pop=sum(-D_pop.*log2(D_pop),'omitnan');
    
 
    FIR=[1 1 1;1 9 1;1 1 1];
    elenum=sum(FIR(:))+1;
%     transi=S_input+(sum(FIR(:))+1)*S_pop;
%     transi=transi(x,y);
%     D_transi=hist(transi(:),(1:2*elenum)-1);

    S_input=conv2(single(cells),repmat(FIR,1,1,size(cells,3)),'same');
    S_input=S_input(x,y);   
    D_input=hist(S_input(:),(1:elenum)-1);
    D_input=D_input/n^2;
    H_input=sum(-D_input.*log2(D_input),'omitnan');
if stepnumber>2    
    Tinput=S_input+S_input_oold.*elenum;
    D_Tinput=hist(Tinput(:),(1:elenum^2)-1);
    SD_Tinput=reshape(D_Tinput,elenum,[]);
    
    VSD_Tinput=SD_Tinput./repmat(sum(SD_Tinput,1),length(SD_Tinput),1);
    VSD_Tinput=-VSD_Tinput.*log2(VSD_Tinput);
    VSD_Tinput(isnan(VSD_Tinput))=0;
    
    HSD_Tinput=SD_Tinput./repmat(sum(SD_Tinput,2),1,length(SD_Tinput));
    HSD_Tinput=-HSD_Tinput.*log2(HSD_Tinput);
    HSD_Tinput(isnan(HSD_Tinput))=0;
    
    D_Tinput=log2(D_Tinput+0.000001);
    
    HD_input=-D_input.*log2(D_input);
    HD_input(isnan(HD_input))=0;
    
    V_input=[HD_input' sum(VSD_Tinput,1)' sum(HSD_Tinput,2) ((1:18))']; 
    [Flux,Udi_Flux]=calcflux(V_input,SD_Tinput);
    V=V_input(:,[2,3 ,1])';
    VV=combvec(V,V);
    VVr=padarray(VV(1:2,:),1,0,'post');
    VVc=padarray(VV(4:5,:),1,0,'post');
    dVV=[VV(1:3,:)-VV(4:6,:)];
    wt=SD_Tinput(:)';
    pdata=[VV(4:6,:);dVV];
    pdata=pdata(:,wt>10);
    NFlux=-cross(VVr,VVc);
    NFlux=NFlux(3,:);
    NFlux=NFlux.*SD_Tinput(:)';
    NFlux=sum(NFlux);
    ANFlux=sum(abs(NFlux));
%     ANFlux=log2(ANFlux);
%     NFlux=log2(abs(NFlux)+1)*sign(NFlux);
    RFlux=sum(VV.^2,1)/2;
    RFlux=RFlux.*SD_Tinput(:)';
    RFlux=sum(RFlux);
%     RFlux=log2(abs(RFlux)+1)*sign(RFlux);

%     WtdFlux=[];
%     LSD_Tinput=SD_Tinput(:);
%     for i=1:length(LSD_Tinput)
%         WtdFlux=[WtdFlux repmat(Flux(:,i),1,LSD_Tinput(i))];
%     end
%     
%     Di_Flux=sum(sum(Flux,3),2)';
% WtdFlux=(VVc-VVr)
%     [theta,rho]=cart2pol(WtdFlux(1,:),WtdFlux(2,:));
%         [T,R]=rose(theta);
%         R=log2(R);
%     Di_Flux=log2(abs(Di_Flux)).*sign(Di_Flux);
%     Abs_Udi_Flux=log(sum(Udi_Flux(:)));
%   UBM_SD_Tinput=D_Tinput(1:elenum);
%     MeVSD_Tinput=VSD_Tinput(1:elenum/2,:)+VSD_Tinput(elenum/2+1:end,:);
%     MeHSD_Tinput=HSD_Tinput(1:elenum/2,:)+VSD_Tinput(elenum/2+1:end,:);

%     LMSD_Tinput=MSD_Tinput(:,1:elenum/2);
%     RMSD_Tinput=MSD_Tinput(:,elenum/2+1:end);
%  
%     MeVSD_Tinput=MeVSD_Tinput(:)';
%   
    wt=SD_Tinput(:)';
    v=[D_input' V_input(:,2)];
    v1=sum(dot(v,v,2));
    wtdVV=dVV(1:2,:).*repmat(wt,2,1);
    v3=mean(wtdVV,2)';
    v2=sum(sum(wtdVV.^2,1).^0.5);
%     v2=sum(sum(v));
%     

%     dmat_V_input=squareform(pdist(V_input(:,1:3)));
%     E_dmat_V_input=[mean(V_input(:,[2 3 1]),1) eigs(dmat_V_input)'];
%     
% H_bool=(V_input(:,2)>1)+2*(V_input(:,3)>1);
% D=zeros(2,2);
% for i=[0 1]
%     for j=[0 1]
%         D(i+1,j+1)=sum(2.^(V_input(H_bool==i+2*j)),1);
%     end
% end
% adj=SD_Tinput;
% A=zeros(4,4);
% for i=0:3;
%     for j=0:3;
%         r=V_input(H_bool==i,4);
%         c=V_input(H_bool==j,4);
%         
%         [r,c]=meshgrid(r,c);
%         cmb=cat(2,r,c);
%         cmb=reshape(cmb,[],2);
% 
%         indices=sub2ind(size(adj), cmb(:,1), cmb(:,2));
% 
%         A(i+1,j+1)=sum(adj(indices));
%     end
% end
% % A=log(A);
% % D=log(D);
% % sum(A(:))
% % sum(D(:))
% L_A=A(:)'/sum(A(:));
% L_D=D(:)'/sum(D(:));
    
    for i=1:length(scalars);
    scalar=scalars{i};
    eval(sprintf('%so=[%so;%s];',scalar,scalar,scalar))
    end
end    

S_input_oold=S_input_old;

S_input_old=S_input;
if update1==1 && stepnumber~=1
    if stepnumber==2
    figure(2)
%     clf
    subplot(2,2,1)
    f1=imagesc(S_input);
    colorbar      
%     legend(cellstr(num2str((1:18)')))
    subplot(2,2,3)
    colormap(gca,'colorcube')
    f3q=quiver(pdata(1,:),pdata(3,:),pdata(4,:),pdata(6,:),0);
%     f3q=quiver3(pdata(1,:),pdata(2,:),pdata(3,:),pdata(4,:),pdata(5,:),pdata(6,:),0);

    markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
    marker=cell2mat(markers([1:9 1:9]));
    f3=gca;
    xlim([0,4])
%     ylim([0,4])
    ylim([0 0.75])
    hold on
    f3g=gscatter(V(1,:),V(2,:),1:18,colormap(gca),marker);
    legend(cellstr(num2str((-1:17)')),'Location','bestoutside')
    hold off


     caxis([0 1000])
    subplot(2,2,4)

    colorbar
    caxis([0,4])
    view(-90,90)
    
    subplot(2,2,2)
    f2=scatter3(V_input(:,2),V_input(:,3) ,V_input(:,1),5,(1:18)*3,'x');
    
    colormap(gca,'colorcube')
    hold on
    axis([0 4 0 4 0 0.5 ])
    view(50,90);

    xlabel('Forward uncertainty')
    ylabel('Backward uncertainty')
    zlabel('Entropic weight')
    rind=floor((k-1)/soupmax)+1;
    suptitle(['MISC' ' rulenumber='  num2str(rind) ' ' rulename{rind,:}])
    elseif mod(stepnumber,2)==0 
        set(f1,'CData',S_input)
        set(f3q,'XData',pdata(1,:),'YData',pdata(3,:),'UData',pdata(4,:),'VData',pdata(6,:));
%         set(f3q,'XData',pdata(1,:),'YData',pdata(2,:),'ZData',pdata(3,:),'UData',pdata(4,:),'VData',pdata(5,:),'WData',pdata(6,:));
    
        for i=1:length(f3g);
            set(f3g(i),'XData',V(1,i),'YData',V(3,i));
%             set(f3g(i),'XData',V(1,i),'YData',V(2,i),'ZData',V_input(i,1));
        end
%           quiver(f3,pdata(1,:),pdata(2,:),pdata(3,:),pdata(4,:),0);
%           hold on
%             f3g=gscatter(V(1,:),V(2,:),1:18,colormap(gca),marker);
%              hold off
%           f3.XLim=[0 4];
%           f3.YLim=[0 4];
%           hold(f3)
%           subplot(f3)
%           ,'filled',    'CData',1:18 
%           hold(f3)
%           scatter(f3,V(1,:),V(2,:),cellstr(num2str((1:18)')))
%           xlim([0,3])
%          ylim([0
% pause    ,4])
%    set(f2,'CData',A);
%         subplot(2,2,3)
%         rose(theta)

%         PP=polarplot(T,R);
%         set(f3,'RData',R,'ThetaData',T);
%         PP.RData
%         PP=gca;
%         ylim([0,50])
%          rlim([0 20])
%          hold on
%         set(f4,'CData',D);
%       set(f3,'CData',dmat_V_input);
        
        subplot(2,2,2)
        scatter3(V_input(:,2),V_input(:,3) ,V_input(:,1),5,(1:18)*3,'x');
                                view(50,90);
                    colormap(gca,'colorcube')
                        view(50,90);
     hold on
%     delete(axs(inds==lag))
%     inds(inds==lag)=[];

%     XYZ=V_input(:,1:3);
%     adj=reshape(D_Tinput,18,[])+20;
%     ax=make3dgraph(adj,XYZ,thres);
%     axs=[axs;ax];
%     inds=[inds+1;ones(length(axs),1)];
%     axis([0 4 0 4 0 0.75 ])
%     colorbar

%      view(50,37.5);
%         set(f4,'XData',[f4.XData V_input(:,2)'],+...
%             'YData',[f4.YData V_input(:,3)'],+...
%             'ZData',[f4.ZData V_input(:,1)'],+...
%             'CData',[f4.CData (1:18)*3]);
    end
    drawnow
%     stem(D_transi)
%     ylim([0,1000])
%     xlim([0 2*(sum(FIR(:))+1)])
end

% fprintf(' active=%d',sum(D_transi~=0))
if pause1==1
pause
end
% pause    
end
% subplot(2,2,2)
% hold off
% 
% subplot(2,2,3)
% hold off

% if store==1
%     dir='gallery/ents/';
%     s=rulename(rind);
%     s=strrep(s,'/','');
%     figname=sprintf('%s_soup%d.jpg',s{1},k);
% 
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
%     saveas(gcf,[dir figname])
% end
inds=[inds;k];

Dv1o=detrend(v1o);
timp_v1o=1-var(Dv1o)/var(v1o);
Mv1o=mean(v1o);
v1s=[v1s;Mv1o timp_v1o];
fprintf('k=%d  Mv1o=%.3f timp_v1o=%.4f',k,Mv1o,timp_v1o)

Dv2o=detrend(v2o);
timp_v2o=1-var(Dv2o')/var(v2o');
Mv2o=mean(v2o);
v2s=[v2s;Mv2o timp_v2o];

% smD_transio=movmean(D_transio,6);
% dmat_transio=squareform(pdist(smD_transio,'euclidean'));

% CD_transio=cov(D_transio);

% D_inputo=movmean(D_inputo,6,'EndPoints','discard');
% DD_inputo=detrend(D_inputo);
% DD_inputo=D_inputo;
% DD_inputo=movmean(DD_inputo,6);
% DD_inputo(1,:)=[];
% timp_D_input=1-var(DD_inputo,1)./var(D_inputo);
% timp_D_inputo=[timp_D_inputo;timp_D_input];
% DH_popo=detrend(H_popo);
% timp_H_pop=1-var(DH_popo,1)./var(H_popo);
% timp_H_popo=[timp_H_popo;timp_H_pop];

% CD_inputo=cov(D_inputo);
% ind=CD_inputo~=0;
% try
% CD_inputo=reshape(CD_inputo(ind),sum(ind(:))^0.5,[]);
% CCD_inputo=corrcov(CD_inputo);
% CCD_inputo=(CCD_inputo+CCD_inputo')/2;
% catch
% CCD_inputo=zeros(size(ind));
% CD_inputo=zeros(size(ind));
% end
% [egv,eg]=eig(CCD_inputo);
% eg=diag(eg);
% eg(eg<0)=0;
% H_eg=sum(-eg.*log2(eg),'omitnan');
% egs=[egs padarray(eg,18-length(eg),0,'pre')];
% eso=[eso mean(CCD_inputo(:))];
% H_ego=[H_ego H_eg];

[pc,score,latent,tsquare] = pca(D_inputo);
% E_dmat_V_inputos=[E_dmat_V_inputos;E_dmat_V_inputo];
if update2==1
makeplot
end
if pause2==1
pause
end
% DH_inputo=detrend(H_inputo);
% DH_inputo=DH_inputo(2:end);
% [p,h]=ksdensity(DH_inputo);
%  plot(h,p)
%  DD_inputo=detrend(D_inputo);
%  p=p/sum(p)
%  H_p=sum(-p.*log2(p));
%  title(H_p)
%  po=[];
%  ho=[];
 
%  for i=1:length(DD_inputo)
% 
%  [p,h]=ksdensity(DD_inputo(:,i));
%  p=p/sum(p)
%  H_p=sum(-p*log2(p));
 
%  end


%  waterfall(DD_inputo')
 
 % [V D] = eigs(CCD_inputo,2, 'SA');
% plot(sort(V(:,2)), '.-')
% ylim([0 1])
% eso=[eso;mean(CCD_inputo(:))];
% imagesc(pinv(CD_inputo))


end
% %%
% 
% % subplot(2,2,4)
% 
% % 
% % M_timp_D_inputo=mean(timp_D_inputo,2,'omitnan');
% m=length(egs);
% xss=[];
% yss=[];
% lss=[];
% for rind=1:floor((m-1)/soupmax)+1; 
%     xs=[];
%     ys=[];
%     for k=((rind-1)*3+1:(rind)*3)
%         if k<=m
%         xs=[xs egs(:,k)];
% %         ys=[ys eso(k)];
%         ys=[ys H_ego(k)];
% %         lss=[]
% %         lss=[lss size(S_poplist{k},1)];
%         end
%     end
%     xss=[xss mean(xs,2)];
%      yss=[yss;mean(ys)];
% end
% % 
% lss=lss/mean(lss);
% subplot(2,2,1)
% plot(eso)
% hold on 
% plot(lss)
% hold off
% title('sum cov (soupwise)')
% subplot(2,2,2)
% plot(timp_H_popo)
% title('variation explained by trend in density-wise entropy (soupwise)')
% plot(yss)
% 
% subplot(2,2,3);
% % scatter(xss(end,:),xss(end-1,:));
% % scatter3(xss(end,:),xss(end-1,:),yss);
% xlabel('1st eigenval');
% ylabel('2nd eigenval');
% title('eig scatter (rulewise)')
% % text(xss(end,:),xss(end-1,:),num2str((1:length(xss))'))
% % text(xss(end,:),xss(end-1,:),yss,num2str((1:length(xss))'))
% 
% subplot(2,2,4);
% imagesc(egs)
% title('heatmap of eigenvalues for each soup')
% colorbar
% 
% suptitle('statistics for rules with similar sum-cov')
% % M_timp_D_inputo=M_timp_D_inputo(~isnan(M_timp_D_inputo)&~isinf(M_timp_D_inputo));
% % plot(M_timp_D_inputo);
% % plot(M_timp_H_pop);
% % hold on
% 
% % plot(M_timp_D_inputo'-detrend(M_timp_D_inputo'))
% % subplot(2,2,4)
% % plot(std(timp_D_inputo','omitnan'))
% hold off

%
% T=table((1:length(rulename))',rulename,yss);
% writetable(T,'rsumcov.xls')