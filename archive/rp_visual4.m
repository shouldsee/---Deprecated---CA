
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
scalars={'D_input','H_input','sH_input','S_inputp',...
    'D_pop_lagged','HD_pop_lagged',+...
    'Hx2','Hxy2','HH_Tinput',...
    'HD_linput','D_linput','MI','sMI','MI_input',...
     'HS_vic','H_pop','D_Tinput','VSD_Tinput','V_input','v3','va','vr','H_Tinput',...
     'LS_pop','LS_input','LS_linput',...
     'KL_div'};
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

if ~exist('lag')
        lag=1;
end
if exist('LS_poplist','var')
    S_poplist=LS_poplist;
end
preallocate


for k=k0:length(S_poplist)
    for i=1:length(scalars)
        eval([scalars{i} 'o=[];'])
       
    end
    S_input_old=zeros(n,n,lag+1);
    S_pop_old=zeros(n,n,lag+1);
    S_linput_old=zeros(n,n,lag+1);
    S_pop_lagged=zeros(n,n);
    S_pop_lagged_old=zeros(n,n);
    S_linput=zeros(n,n);
    SD_Tlinput=zeros(n,n);
    
    H_Tinput=zeros(1,lag+1);
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
    LS_pop=S_pop(:)';

    cells(x,y)=S_pop;
    cells(1,:)=cells(x(end),:);
    cells(end,:)=cells(x(1),:);
    cells(:,1)=cells(:,y(end));
    cells(:,end)=cells(:,y(1));
    
    D_pop=hist(S_pop(:),0:1);
    D_pop=D_pop/sum(D_pop);
    H_pop=sum(-D_pop.*log2(D_pop),'omitnan');
    
    
%     FIR=[1 1 1;1 9 1;1 1 1];
%     elenum.S_input=sum(FIR(:))+1;
%     transi=S_input+(sum(FIR(:))+1)*S_pop;
%     transi=transi(x,y);
%     D_transi=hist(transi(:),(1:2*elenum)-1);

    S_input=conv2(single(cells),repmat(FIR.S_input,1,1,size(cells,3)),'same');
    S_input=S_input(x,y);
%     [CS_input,base]=autoMI(S_input);
%     CS_input=CS_input/base;
%     LCS_input=CS_input(:)';
    LS_input=S_input(:)';
    S_input_old=cat(3,S_input,S_input_old(:,:,1:lag));
    
%     S_inputo=reshape(LS_inputo,size(LS_inputo,1),30,30);
    D_vic=zeros(1,elenum.S_input^2);
%     D_vic=[];
    for i=1:numel(FIR.vic)
        D_vic=D_vic+hist(reshape(convn(S_input_old,FIR.vic{i},'valid'),1,[]),(1:elenum.S_input^2)-1);
    end
    
    SD_vic=reshape(D_vic,numel(D_vic)^0.5,[]);
%     SD_vic=(SD_vic+SD_vic');
    HSD_vic=entropise(SD_vic);
    HS_vic=sum(HSD_vic(:));
    try
    LS_input_old=[LS_input;LS_inputo(1:lag,:)];
    catch
    LS_input_old=[LS_input;zeros(lag,n^2)];
    end
    
    input=S_input;
    CS_inputX=circenum(input,[1 0]);
    CS_inputY=circenum(input,[0 1]);
    CS_inputX=reshape(CS_inputX,size(CS_inputX,1),[]);
    CS_inputY=reshape(CS_inputY,size(CS_inputY,1),[]);
    
    Dxy2=sum(sum(vcalc_Dxy(CS_inputX([1:2 30],:),CS_inputY(1:2,:)),1),2);
    SDxy2=reshape(Dxy2,numel(Dxy2)^0.5,[]);
    HSDxy2=entropise(SDxy2);
    Dx2=sum(SDxy2,1);
    HDx2=entropise(Dx2);
    Hx2=sum(HDx2(:));
    HDxy2=entropise(Dxy2(:));
    
    Hxy2=sum(HDxy2);

%     dCS_input=pdist2(CS_inputX,CS_inputY,@calc_MI);
%     dCS_input=circshift(full(dCS_input),[15 15]);
%     dCS_input=dCS_input/dCS_input(15,15);
%     ldCS_input=real(log2(dCS_input));

    
    
    D_input=hist(S_input(:),(1:elenum.S_input)-1);
    D_input=D_input/sum(D_input(:));
    bn=binopdf(0:8,8,D_pop(2));
    M_input=[bn*D_pop(1) bn*D_pop(2)];
    KL_div=sum(D_input.*log2(D_input./M_input),'omitnan');
    S_input_V=permute(18.^(0:size(S_input_old,3)-1),[1,3,2]);
    S_input_otemp=bsxfun(@times,S_input_old,S_input_V);
    sHD_input=D_input.*log2(D_input).^2;
    sHD_input(isnan(sHD_input))=0;
    sH_input=sum(sHD_input);
    HD_input=-D_input.*log2(D_input);
    HD_input(isnan(HD_input))=0;
    H_input1=sum(HD_input(:),'omitnan');
    sH_input=sH_input-H_input1^2;
    
%     D_input2=hist(reshape(sum(S_input_otemp(:,:,1:2),3),1,[]),0:sum(17*S_input_V(:,:,1:2)));
%     D_input2=D_input2/sum(D_input2(:));
%     HD_input2=-D_input2.*log2(D_input2);
%     H_input2=sum(HD_input2(:),'omitnan');
try
    D_input3=hist(reshape(sum(S_input_otemp(:,:,1:3),3),1,[]),0:sum(17*S_input_V(:,:,1:3)));
    D_input3=D_input3/sum(D_input3(:));
catch
    D_input3=ones(1,n);
end
    HD_input3=-D_input3.*log2(D_input3);

    H_input3=sum(HD_input3(:),'omitnan');
    
%     H_input=[H_input1 H_input2 H_input3];
    H_input=H_input1;
    S_pop_old=cat(3,S_pop,S_pop_old(:,:,1:lag));
%     S_linput_old=cat(3,S_linput,S_linput_old(:,:,1:lag));
%     elenum=elenum;
    S_pop_lagged_old=S_pop_lagged;
    S_pop_lagged=sum(S_pop_old.*repmat(pop_elenumV,n,n,1),3);
    D_pop_lagged=hist(reshape(S_pop_lagged(1:30,2:29),1,[]),(1:elenum.S_pop_lagged)-1);
    D_pop_lagged=D_pop_lagged./sum(D_pop_lagged);
    HD_pop_lagged=entropise(D_pop_lagged);
    H_pop_lagged=sum(HD_pop_lagged);
%     FIR=reshape(elenum.S_pop_lagged.^([0 1 2 3]),2,2);
%     FIR=reshape(elenum.S_pop_lagged.^([0]),1,1);
%     elenum.S_linput=sum(FIR(:).*(elenum.S_pop_lagged-1))+1;
%     FIR.S_linput=elenum.S_pop_lagged.^[0 1;2 3];
% %     FIR.S_linput=elenum.S_pop_lagged.^[0 1;2 3];
%     S_linput=conv2(S_pop_lagged,FIR.S_linput,'same');
%     LS_linput=S_linput(:)';
%     %     S_linput=dict.S_linput(1+S_linput);
% %     elenum.S_linput=max(dict.S_linput);
% %     elenum.S_linput=sum(FIR.S_linput(:).*(elenum.S_pop_lagged-1))+1;     
% %     S_Tlinput=S_linput+S_linput_old(:,:,la+1).*elenum.S_linput;
% %     S_Tinput_lagged=S_input_lagged
%     D_linput=hist(S_linput(:),(1:elenum.S_linput)-1);
%     D_linput=D_linput/sum(D_linput);
%     HD_linput=entropise(D_linput);
%     H_linput=sum(HD_linput(:));
    sMI(1)=(2*H_pop_lagged-H_linput)/H_pop_lagged;


    
    
%     FIR.S_linput=reshape([elenum.S_pop_lagged.^[0 1 2] 0],2,2);
%     S_linput=conv2(S_pop_lagged,FIR.S_linput,'same');
%     elenum.S_linput=sum(FIR.S_linput(:).*(elenum.S_pop_lagged-1))+1;     
%     D_linput=hist(S_linput(:),(1:elenum.S_linput)-1);
%     D_linput=D_linput/sum(D_linput);
%     HD_linput=entropise(D_linput);
%     H_linput2=sum(HD_linput(:));
%     sMI(2)=()
    
    Tinput=S_input+S_input_old(:,:,lag+1).*elenum.S_input;
%     Tinput=S_pop_lagged+pop_elenum.*S_pop_lagged_old;
%     D_Tlinput=hist(Tlinput(:),(1:elenum.S_linput^2)-1);
%     SD_Tlinput=reshape(D_Tlinput,elenum.S_linput,[])';
%     SD_Tinput=SD_Tlinput;
    D_Tinput=hist(Tinput(:),(1:elenum.S_input^2)-1);
%     SD_Tinput=reshape(D_Tinput,elenum.S_input,[])';
    D_Tinput=D_Tinput/sum(D_Tinput);
    SD_Tinput=reshape(D_Tinput,elenum.S_input,[])';
    HHSD_Tinput=entropise(SD_Tinput);
    HH_Tinput=sum(HHSD_Tinput(:));
    SD_Tinput_marg1=sum(SD_Tinput,1);
    SD_Tinput_marg2=sum(SD_Tinput,2);
    
    HD_Tinput=-D_Tinput.*log2(D_Tinput);
    HD_Tinput(isnan(HD_Tinput))=0;
    H_Tinput=sum(HD_Tinput);
    
    
    sHV=sum(entropise(sum(SD_Tinput,1)));
    sHH=sum(entropise(sum(SD_Tinput,2)));

    VSD_Tinput=SD_Tinput./repmat(sum(SD_Tinput,1),length(SD_Tinput),1);
    VSD_Tinput=-VSD_Tinput.*log2(VSD_Tinput);
    VSD_Tinput(isnan(VSD_Tinput))=0;
    HSD_Tinput=SD_Tinput./repmat(sum(SD_Tinput,2),1,length(SD_Tinput));
    HSD_Tinput=-HSD_Tinput.*log2(HSD_Tinput);
    HSD_Tinput(isnan(HSD_Tinput))=0;    
    HV=sum(VSD_Tinput,1)';
    HH=sum(HSD_Tinput,2);
    MI_input=sHV+sHH-H_Tinput;
    
    
    elenum.plot=getfield(elenum,'S_input');
if stepnumber>lag+1 % && mod(stepnumber-1,lag)==0
    
    Tlinput=S_linput+S_linput_old(:,:,lag+1).*elenum.S_linput;
    la=lag;
    SD_linput=SD_Tinput;
    MI=-SD_Tinput.*log2(bsxfun(@times,SD_Tinput_marg1,SD_Tinput_marg2)./SD_Tinput);
    MI=sum(MI(:),'omitnan');
   
      D_Tinput=log2(D_Tinput+0.000001);
    
      HD_input=zeros(1,length(HH));
    V_input=[HD_input' HV+HH HV-HH ((1:elenum.plot))']; 
%      V_input=[HD_pop_lagged' HV+HH HV-HH]; 

    [Flux,Udi_Flux]=calcflux(V_input,SD_Tinput);
    V=V_input(:,[2,3 ,1])';
    VV=combvec(V,V);
    VVr=padarray(VV(1:2,:),1,0,'post');
    VVc=padarray(VV(4:5,:),1,0,'post');
    dVV=[VV(4:6,:)-VV(1:3,:)];
    wt=SD_Tinput(:)';
    wt=wt/sum(wt);
%     wt=-wt.*log2(wt);
%     wt(isnan(wt))=0;
    pdata=[VV(1:3,:);dVV];
    pdata=pdata(:,wt>0.01);  
%     wt=SD_Tinput(:)';
%     v=[D_input' V_input(:,2)];
%     v1=sum(dot(v,v,2));
    wtdVV=dVV(1:2,:).*repmat(wt,2,1);
    v2=mean(wtdVV,2)';
    va=mean(abs(wtdVV),2)';
    vr=mean(wtdVV,2)';
    
    S_inputp=permute(S_input,[3 1 2]);
    for i=1:length(scalars);
        scalar=scalars{i};
        eval(sprintf('%so=[%so;%s];',scalar,scalar,scalar))
    end
end    

% S_input_oold=S_input_old;

% S_input_old=S_input;
if update1==1 && stepnumber>lag+1
    if stepnumber==lag+2
    figure(2)
%     clf
    subplot(2,2,1)
    f1=imagesc(S_linput);
    colorbar      
    caxis([0 elenum.S_input])
    view([90 -90])
%     legend(cellstr(num2str((1:elenum)')))
    subplot(2,2,3)
    colormap(gca,'colorcube')
    f3q=quiver(pdata(1,:),pdata(3,:),pdata(4,:),pdata(6,:),0);
%     f3q=quiver3(pdata(1,:),pdata(2,:),pdata(3,:),pdata(4,:),pdata(5,:),pdata(6,:),0);

    markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
    marker=cell2mat(markers([1:9 1:9]));
    f3=gca;
    xlim([0,4])
%     ylim([0 0.75])
    ylim([-4 4])
    zlim([0 0.75])
    
    hold on
    f3g=gscatter(V(1,:),V(2,:),1:elenum.plot,colormap(gca),marker);
    legend(cellstr(num2str((0:elenum.plot)')),'Location','bestoutside')
    hold off
    refline([0 0])
    caxis([0 1000])
    subplot(2,2,4)
%     f4i=imagesc(VSD_Tinput);
%     f4i=imagesc(ldCS_input);
    f4i=imagesc(HSD_vic);
    colorbar
%     caxis([0,4])
    view(-90,90)
%     ax4=gca;
%     ax4.XTickLabels=num2str((ax4.XTick-1)')';
    
    subplot(2,2,2)
    f2i=imagesc(HHSD_Tinput);
    view(-90,90)
%     f2=scatter3(V_input(:,2),V_input(:,3) ,V_input(:,1),5,(1:elenum.plot)*3,'x');
%     
%     colormap(gca,'colorcube')
%     hold on
%     axis([0 4 0 4 0 0.5 ]);
%     view(50,90);
% 
%     xlabel('Forward uncertainty')
%     ylabel('Backward uncertainty')
%     zlabel('Entropic weight')
    rind=floor((k-1)/soupmax)+1;
    suptitle(['MISC' ' rulenumber='  num2str(rind) ' ' rulename{rind,:}])
    elseif mod(stepnumber,lag)==0 
        set(f1,'CData',S_input)
        set(f2i,'CData',HHSD_Tinput)
%         set(f1,'CData',S_linput);
%         set(f3q,'XData',pdata(1,:),'YData',pdata(3,:),'UData',pdata(4,:),'VData',pdata(6,:));
        set(f3q,'XData',pdata(1,:),'YData',pdata(2,:),'ZData',pdata(3,:),'UData',pdata(4,:),'VData',pdata(5,:),'WData',pdata(6,:));
        set(f4i,'CData',HSD_vic);
        for i=1:length(f3g);
%             set(f3g(i),'XData',V(1,i),'YData',V(3,i));
            set(f3g(i),'XData',V(1,i),'YData',V(2,i),'ZData',V_input(i,1));
        end
%         subplot(2,2,2)
%         scatter3(V_input(:,2),V_input(:,3) ,V_input(:,1),5,(1:elenum.plot)*3,'x');
% %                                 view(50,90);
%                     colormap(gca,'colorcube')
%                            view(50,90);
%     hold on
    end
    drawnow
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
% 
% Dv1o=detrend(v1o);
% timp_v1o=1-var(Dv1o)/var(v1o);
% Mv1o=mean(v1o);
% v1s=[v1s;Mv1o timp_v1o];
% fprintf('k=%d  Mv1o=%.3f timp_v1o=%.4f',k,Mv1o,timp_v1o)
% 
% Dv2o=detrend(v2o);
% timp_v2o=1-var(Dv2o')/var(v2o');
% Mv2o=mean(v2o);
% v2s=[v2s;Mv2o timp_v2o];


[pc,score,latent,tsquare] = pca(D_inputo);

if update2==1 && stepnumber>20
makeplot
end
if pause2==1
pause
end
end