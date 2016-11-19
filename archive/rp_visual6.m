pre_init
preallocate
% ruletemp=importrule(rulename);
kmax=length(rulename)*soupmax;
if randrule
    kmax=100000000;
end

for k=k0:kstep:kmax

    initialise
    
if sum(S_popo(:))==0
        for i=1:length(lists)
           list=lists{i};
           eval(sprintf('%ss=[%ss;zeros(1,size(%ss,2))];',list,list,list))
        end
    
        continue
end

rind=floor((k-1)/soupmax)+1;    
    
    if ~randrule
    r=ruletemp(rind,:);
    rule=[r{1},r{2}];
    else
    ktemp=k;
    k=8;
    rrule=1;
    update_rule;
    ruletemp(rind,:)={ruleB,ruleS};
    rulename{rind}=rstr;
    k=ktemp;
    rule=rulecurr;
%     rstr
%     pause
    end
    rulecurr=rule;
%     B,S=
%     [ruleB,ruleS]=BSrule(rule)

for stepnumber=1:stepmax    S_pop=cells(x,y);

%     size(S_popo,1)    
    fprintf('\n clip k=%d,step %d ',k,stepnumber)
%     S_pop=squeeze(S_popo(stepnumber,:,:));
%     S_pop=cells(x,y);
    S_pop_old=S_pop;
    LS_pop=S_pop(:)';
    S_pop=cells(x,y);
    
if stepnumber>1
end
    if ~allrand
        cells(x,y)=rulecurr(S_input+1);
    else
%         mx=mod(stepnumber,2*length(po));
%         md=floor(mx/2);
        md=mod(stepnumber,length(po));
%         if md==mx/2;
%         pcurr=po(10+1);
        pcurr=po(md+1);
%            cells(x,y)=rand(length(x),length(y))<;
          cells(x,y)=rand(length(x),length(y))<pcurr;
          S_pop=cells(x,y);
          cellsM=repmat(cells,[1,1,uni]);
          mt=randi([1,n^2],[uni,1]);
          [a,b]=ind2sub(size(cells),mt);
          cinds=sub2ind(size(cellsM),a',b',1:size(cellsM,3));
          cellsM(cinds)=1-cellsM(cinds);
%           'respawn'
%         else    
%         end

       for tiny=1:warmup; 
            S_inputT=conv2(cells,FIR.S_input,'same');
            S_input=S_inputT(x,y);
            cells(x,y)=rulecurr(S_input+1);                                    
            cells(1,:)=cells(x(end),:);
            cells(end,:)=cells(x(1),:);
            cells(:,1)=cells(:,y(end));
            cells(:,end)=cells(:,y(1));
            if exist('cellsM','var')
%                 cellsM(x,y,:)=S_pop;
                cellsM(1,:,:)=cellsM(x(end),:,:);
                cellsM(end,:,:)=cellsM(x(1),:,:);
                cellsM(:,1,:)=cellsM(:,y(end),:);
                cellsM(:,end,:)=cellsM(:,y(1),:);
                S_inputMT=convn(single(cellsM),FIR.S_input,'same');
                S_inputM=S_inputMT(x,y,:);
                cellsM=rulecurr(S_inputMT+1);
                LS_popM=squeeze(reshape(cellsM(x,y,:),n^2,1,[]))';
            end
            LS_pop=reshape(cells(x,y),1,[]);;
%             d=pdist2(LS_pop,LS_popM,'hamm');
            
            d=mean(bsxfun(@ne,LS_pop,LS_popM),2)';
%             set(h,'Data',d);
            drawnow
%             pause(0.2)
            LS_input_log(tiny,:)=S_input(:)';
            dall(tiny,:)=d;
            if update1
                makeupdate
            end
       end

    end
cells(1,:)=cells(x(end),:);
cells(end,:)=cells(x(1),:);
cells(:,1)=cells(:,y(end));
cells(:,end)=cells(:,y(1));

    D_pop=hist(S_pop(:),0:1);
    D_pop=D_pop/sum(D_pop);
    H_pop=sum(-D_pop.*log2(D_pop),'omitnan');
    
%     S_inputPT=conv2(single(xor(cells,S_defectT)),FIR.S_input,'same');
%     S_inputP=S_inputPT(x,y);
%     
    S_inputT=conv2(single(cells),repmat(FIR.S_input,1,1,size(cells,3)),'same');
    S_input=S_inputT(x,y);

%     D_input=hist(S_input(:),(1:elenum.S_input)-1);
%     D_input=D_input/sum(D_input(:));
    
    LS_input=S_input(:)';
    LS_input_old=cat(1,LS_input,LS_input_old(1:lag,:));
    
%     St_input_old=bsxfun(@times,elenum.S_input.^(0:lag-1)',LS_input_old(1:lag,:));
%     S_input_old=cat(3,S_input,S_input_old(:,:,1:lag));
%     D_input_old=shiftdim(sum(reshape(DS_input_old,[],elenum.S_input),1),1);
%     H_input_old=sum(entropise(D_input_old(:)));   

%     D_input_old=hist(S_input_old(:),(1:elenum.S_input)-1);
%     P_input_old=D_input_old/sum(D_input_old(:));

%     ID_input=(LS_input~=LS_input_old(2,:));
%     ID_input=find(ID_input);
%     ID_input=1:length(LS_input);
%     if ~isempty(ID_input) 
%     MI_input=pdist(LS_input_old,@calc_MI);
% %     MI_input=pdist(LS_input_old,'hamming');
% %     calc_MI(LS_input(ID_input),LS_input_old(2,ID_input));
%     else
%     MI_input=0;
%     end

%     
%     MI_logV=calc_MI(LS_input_log(1,:),LS_input_log)';
%     H_logall=calc_Hy(LS_input,LS_input_log);
%     H_logall=pdist2(LS_input,LS_input_log,'hamm');
%   
%     for i=1:size(LS_input_log,1)
%        H_logall(i)=size(dzip(reshape(LS_input_log(i,:),n,[])),1);
%     end
%       MI_logall=(squareform(pdist(LS_input_log,@calc_MI)));
%     MI_logall=MI_logall+diag(H_logall);
%     [gd1,gd2]=ndgrid(1:size(MI_logall,1),1:size(MI_logall,2));

%     diag(gd1)=10:-1:1
%     MI_input=-log2(pdist(LS_input_old,'hamming'));
%     MI_input(isinf(MI_input))=0;
%     rec=[md MI_input];
%     rec=cat(3,repmat(po(md+1),1,warmup),1:warmup,MI_logV);
%    rec=cat(3,gd1(:)',gd2(:)',MI_logall(:)');
    wv=repmat((1:warmup)',1,uni);
    pv=repmat(pcurr,warmup,uni);
%     rec=cat(3,repmat(pcurr,1,warmup,uni),1:warmup,H_logall(:)');
    rec=cat(3,pv,wv,dall);
%     set(f,'CData',dall)
   %     prec=repmat(po(md+1),1,size(rec,2));
    % make_info
%     make_vic
%     make_flow
%     HLE=std(lossexponent(:));

%     S_inputV=bsxfun(S_input)
%     repS_input=repmat(S_input,1)
%     N_defect=2*sum(S_defect(:)>0)/sum(reshape(S_pop_old(:,:,1:2),1,[]));
%     N_defect=sum(S_defect(:)>0)/900;
    
%     P_input=2^(-H_input);
    
%     LS_defectL=log2(S_defect(:)');
%     S_popND(:,:,1)
%     bsxfun(S_popND,permute(2.^(0:7)',[1 3 2]));
%     S_change=rule_change(S_popND);
%     t1=[t1;toc];    
%       max(S_input(:))
%     change=rule(S_inputP+1)~=rule(S_input+1);
%     S_defect=change;
%     S_defectT(x,y)=S_defect;
%     S_defectT(1,:)=S_defectT(x(end),:);
%     S_defectT(end,:)=S_defectT(x(1),:);
%     S_defectT(:,1)=S_defectT(:,y(end));
%     S_defectT(:,end)=S_defectT(:,y(1));


%     S_pop2input=S_pop+S_input.*elenum.S_pop;
%     D_pop2input=hist(S_pop2input(:),0:elenum.S_pop2input-1);
%     HD_pop2input=entropise(HD_pop2input)
    
    
%     squeeze(reshape(stack(LS_input_old),[],1,n^2));
%     St_input_MI=calc_MI(LS_input,St_input_old);
%     St_input_MI=reshape(St_input_MI,size(St_input_MI,1).^0.5,[]);
%     SV_pop=bsxfun(@eq,S_pop,permute(0:1,[1 3 2]));
%     HD_popx=entropise(sum(SV_pop,1));
%     HD_popy=entropise(sum(SV_pop,2));
%     HD_popxy=entropise(SV_pop);
    
%     MI_pop=sum(HD_popx(:))+sum(HD_popy(:))-sum(HD_popxy(:));
%     S_dinputT=conv2(single(cells),FIR.S_dinput,'same');
%     S_dinput=S_dinputT(x,y);
%     D_dinput=hist(S_dinput(:),0:elenum.S_dinput-1);
%     D_dinput=D_dinput/sum(D_dinput);
%     HD_dinput=entropise(D_dinput);
%     H_dinput=sum(HD_dinput);
%     I_dinput=-log2(D_dinput);
%     I_dinput(isnan(I_dinput))=0;
%     SI_dinputT=I_dinput(S_dinputT+1);
%     SI_dinput=SI_dinputT(x,y);
%     SI_dinput_old=cat(3,SI_dinput,SI_dinput_old(:,:,1:lag));

%     S_Iinput=conv2(SI_dinput,FIR.SI_dinput,'same');
%     LS_dinput=S_dinput(:)';
%     LS_dinput_old=[LS_dinput;LS_dinput_old(1:lag,:)];
%     MIs(1)=H_input;
%     for i=1:3
%         MIs(i+1)=mutualinfo(LS_input_old(1:i+1,:)',1/n^2.*ones(n^2,1));
%     end
%     D_dinput=hist(S_dinput(:),0:elenum.S_dinput-1);
%     HD_dinput=entropise(D_dinput);
%     H_dinput=+sum(HD_dinput(:),1);
%     S_input2dinput=S_input.*1+S_dinput.*elenum.S_input;
    
%     MI_dinput=calc_MI(reshape(LS_dinput_old(1:lag,:),1,[]),reshape(LS_dinput_old(2:lag+1,:),1,[]));

%     D_input2dinput=hist(S_input2dinput(:),0:elenum.S_tot2d-1);
%     D_input2dinput=sum(bsxfun(@eq,S_input2dinput(:),0:elenum.S_tot2d-1),1);
%     HD_input2dinput=entropise(D_input2dinput);
%     H_input2dinput=sum(HD_input2dinput(:),1);
%     H_input2dinput=H_input2dinput-calc_Hxy(S_input(:)',S_dinput(:)');
%     H_input2dinput=+calc_Hxy(S_input(:)',S_dinput(:)');
%     H_input=calc_Hx(S_input(:)',S_dinput(:)');
%     H_dinput=calc_Hy(S_input(:)',S_dinput(:)');
%     H_dinput_input=(H_input2dinput-H_input);




%     make_defect;

% mean(lossexponent(:).^2-mean(lossexpon));
% hedge=0:0.01:2;

% S_flow=ones(n,n);
% S_inputND=findneighbor(S_input,neighbor.S_pop);
% S_flowND=findneighbor(S_flow,neighbor.S_pop);
% % nSMI_vic=SMI_vic/MI_vic;
% bonddict=logsig(SMI_vic);
% id=sub2ind(size(bonddict),S_inputND+1,repmat(S_input+1,1,1,8));
% bond=bonddict(id);
% bond=bsxfun(@rdivide,bond,sum(bond,3));
% bond=makezero(bond);
%     SMI_vic=bsxfun(SMI_vic,)
%     SI_vic=-log2(SP_vic);
%     SI_vic(isnan(SI_vic)|isinf(SI_vic))=0;
%     SI_vic()=0;
%     I_vic=sum((SI_vic.*SP_vic),2)./sum(SP_vic,2);

%     dLSI=-diff(LSI_input_old);
% % f3p=[LSI_input_old(2,:)',dLSI(1,:)'];
% % f3p=[I_input',(MI_vic1)];
% f3p=[SSMI_vic1(:),SSMI_vic1_input(:)];
% f3p=abs(f3p);
% sf3p=sum(f3p,1)/n^2;
% % sf3p(1)-H_dinput;
% FH_input=sf3p(2);
% 


%     SSD_vic=shiftdim(sum(SD_vic,1),1)/2;
%     SSSD_vic=(SSD_vic+SSD_vic')/2;
%     HSD_vic=entropise(SSD_vic);
%     HSSD_vic2=entropise(SSSD_vic);
%     HSS_vic2=sum(HSSD_vic2(:));
%     HS_vic2=sum(HSD_vic(:));
    
%      fHSD_vic=HSD_vic;
  
%     input=S_input;
%     CS_inputX=circenum(input,[1 0]);
%     CS_inputY=circenum(input,[0 1]);
%     CS_inputX=reshape(CS_inputX,size(CS_inputX,1),[]);
%     CS_inputY=reshape(CS_inputY,size(CS_inputY,1),[]);
%     calc_Hy_x=@(a,b)calc_Hxy(a,b)-calc_Hx(a,b);
%     Hy_x=pdist2(CS_inputX,CS_inputY,@calc_MI);
%     Hy_x=Hy_x./Hy_x(1);
    
    
%     x1=Hy_x(2,1).^((1:size(Hy_x,1))-1);
%     x1=x1';
%     y1=Hy_x(1,2).^((1:size(Hy_x,2))-1);
%     GHy_x=bsxfun(@times,x1,y1);
%     Hy_x=(GHy_x-Hy_x)./Hy_x;

%     SDxy2=reshape(Dxy2,numel(Dxy2)^0.5,[]);
%     HSDxy2=entropise(SDxy2);
%     Dx2=sum(SDxy2,1);
%     HDx2=entropise(Dx2);
%     Hx2=sum(HDx2(:));
%     HDxy2=entropise(Dxy2(:));
%     
%     Hxy2=sum(HDxy2);

%     dCS_input=pdist2(CS_inputX,CS_inputY,@calc_MI);
%     dCS_input=circshift(full(dCS_input),[15 15]);
%     dCS_input=dCS_input/dCS_input(15,15);
%     ldCS_input=real(log2(dCS_input));

    
%     
% 
%     bn=binopdf(0:8,8,D_pop(2));
%     M_input=[bn*D_pop(1) bn*D_pop(2)];
%     KL_div=sum(D_input.*log2(D_input./M_input),'omitnan');
%     S_input_V=permute(18.^(0:size(S_input_old,3)-1),[1,3,2]);
%     S_input_otemp=bsxfun(@times,S_input_old,S_input_V);
%     sHD_input=D_input.*log2(D_input).^2;
%     sHD_input(isnan(sHD_input))=0;
%     sH_input=sum(sHD_input);
%     H_input1=sum(HD_input(:),'omitnan');
%     sH_input=sH_input-H_input1^2;
    
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
%     H_input=H_input1;
%     S_pop_old=cat(3,S_pop,S_pop_old(:,:,1:lag));
    S_linput_old=cat(3,S_linput,S_linput_old(:,:,1:lag));
%     elenum=elenum;
    S_pop_lagged_old=S_pop_lagged;
%     S_pop_lagged=sum(S_pop_old.*repmat(pop_elenumV,n,n,1),3);
    D_pop_lagged=hist(reshape(S_pop_lagged(1:30,2:29),1,[]),(1:elenum.S_pop_lagged)-1);
    D_pop_lagged=D_pop_lagged./sum(D_pop_lagged);
    HD_pop_lagged=entropise(D_pop_lagged);
    H_pop_lagged=sum(HD_pop_lagged);
%     sMI(1)=(2*H_pop_lagged-H_linput)/H_pop_lagged;
% 
    Tinput=S_input+S_input_old(:,:,lag+1).*elenum.S_input;
    D_Tinput=hist(Tinput(:),(1:elenum.S_input^2)-1);
    
    P_Tinput=D_Tinput/sum(D_Tinput);
    
    
%     P_Tinput=
%     SD_Tinput=reshape(D_Tinput,elenum.S_input,[]);
%     SP_Tinput=reshape(P_Tinput,elenum.S_input,[]);
%     SI_Tinput=-log2(SP_Tinput);
%     SI_Tinput(isnan(SI_Tinput))=0;
%     SI_Tinput1=sum((SI_Tinput.*SP_Tinput),2)/sum(SP_Tinput,2);
%     HHSD_Tinput=entropise(SD_Tinput);
%     HH_Tinput=sum(HHSD_Tinput(:));
% 
%     HD_Tinput=entropise(D_Tinput);
%     H_Tinput=sum(HD_Tinput);
%     
%     
%     sHV=sum(entropise(sum(SD_Tinput,1)));
%     sHH=sum(entropise(sum(SD_Tinput,2)));
% 
%     MI_input=sHV+sHH-H_Tinput;
%     
    
    elenum.plot=getfield(elenum,'S_input');
if stepnumber>lag+1 % && mod(stepnumber-1,lag)==0
    
%     Tlinput=S_linput+S_linput_old(:,:,lag+1).*elenum.S_linput;
    la=lag;
%     SD_linput=SD_Tinput;
%     MI=-SD_Tinput.*log2(bsxfun(@times,SD_Tinput_marg1,SD_Tinput_marg2)./SD_Tinput);
%     MI=sum(MI(:),'omitnan');
   
%       D_Tinput=log2(D_Tinput+0.000001);
    
%       HD_input=zeros(1,length(HH));
%     V_input=[HD_input' HV+HH HV-HH ((1:elenum.plot))']; 
%      V_input=[HD_pop_lagged' HV+HH HV-HH]; 

%     [Flux,Udi_Flux]=calcflux(V_input,SD_Tinput);
%     V=V_input(:,[2,3 ,1])';
%     VV=combvec(V,V);
%     VVr=padarray(VV(1:2,:),1,0,'post');
%     VVc=padarray(VV(4:5,:),1,0,'post');
%     dVV=[VV(4:6,:)-VV(1:3,:)];
%     wt=SD_Tinput(:)';
%     wt=wt/sum(wt);
%     wt=-wt.*log2(wt);
%     wt(isnan(wt))=0;
%     pdata=[VV(1:3,:);dVV];
%     pdata=pdata(:,wt>0.01);  
%     wt=SD_Tinput(:)';
%     v=[D_input' V_input(:,2)];
%     v1=sum(dot(v,v,2));
%     wtdVV=dVV(1:2,:).*repmat(wt,2,1);
%     v2=mean(wtdVV,2)';
%     va=mean(abs(wtdVV),2)';
%     vr=mean(wtdVV,2)';
    S_inputp=permute(S_input,[3 1 2]);
    for i=1:length(scalars);
        scalar=scalars{i};
        eval(sprintf('%so=[%so;%s];',scalar,scalar,scalar))
    end
end    

% S_input_oold=S_input_old;

% S_input_old=S_input;
% LSI_dinput_old=squeeze(reshape(SI_dinput_old,[],1,lag+1))';
% dLSI=-diff(LSI_dinput_old);
% f3p=[LSI_dinput_old(2,:)',dLSI(1,:)'];
% f3p=abs(f3p);
% sf3p=sum(f3p,1)/n^2;
% sf3p(1)-H_dinput;
% FH_dinput=sf3p(2);

if update1==1 && stepnumber>lag+1
makeupdate
end

% fprintf(' active=%d',sum(D_transi~=0))
if pause1==1
pause
end
% pause    
end

    rind=floor((k-1)/soupmax)+1;

% rind=ceil((k)/soupmax)    
% if rind>size(rulename,1);
%     rind=size(rulename,1)
%     rind=1;
% end
if work==1
stats=[H_inputo FH_inputo];
so.stats=[so.stats;stats];
so.rulelist=[so.rulelist;repmat(rulename(rind,:),size(stats,1),1)];
so.label=[so.label;label*ones(size(stats,1),1)];
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
% inds=[inds;k];
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


% [pc,score,latent,tsquare] = pca(D_inputo);

if update2==1 && stepnumber>20
makeplot
end
if pause2==1
pause
end
end