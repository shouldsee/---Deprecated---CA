window=6;
thres1=0.3;
thres2=0.001;

disp(['effective step=' num2str(estep)]);
if estep>20;
% for i=1:length(scalars)
%     scalar=scalars{i};
% %     fprintf('%d',length(H_inputo))
%     eval(sprintf('%so=%so(30:end,:);',scalar,scalar))
% end

% 
% MD_inputo=movmean(D_inputo,6,'EndPoints','discard');
% % CD_inputo=cov(MD_inputo);
% [coef,~,latent,~,D_PCAV,~]=pca(MD_inputo);
% D_PCAV=D_PCAV/100;
% D_PCAV(D_PCAV<0)=0;
% H_PCAV=sum(-D_PCAV.*log2(D_PCAV),'omitnan');
if exist('MIo','var');
%    out=[H_pop_laggedo,MIo,H_inputo];
%     temp(soupnumber,:)=[mean(out),std(out),estep];
end


if exist('phaseXo','var')
    
% H_linputo=sum(HD_linputo,2);
phaseo=[phaseXo,phaseYo,phaseZo];

mxs=movmean(phaseXo,6,'Endpoints','discard');
mys=movmean(phaseYo,6,'Endpoints','discard');
mzs=movmean(phaseZo,6,'Endpoints','discard');
[coeff,~,latent,~,explained]=pca([mxs,mys]);
wcoeff=bsxfun(@times,latent,coeff')';

% tc=cov(phaseo);
mp=mean([phaseXo,phaseYo],1);
% c=cov(diff(phaseo));
% tv=trace(c);

% sd=std(phaseo,1);
% try
% GMModel = fitgmdist(phaseo,1);
% AIC=GMModel.AIC;
% catch Err
% AIC=nan;
% end
temp(soupnumber,:)=[mp, wcoeff(:)',estep];
% temp(soupnumber,:)=[mp,c(1,1), c(2,1), c(2,2),tv,estep];
end
% if exist('V_inputo','var')
% id=find(H_bool==0);
% V=reshape(V_inputo,[],18,4);
% V=squeeze(mean(V,1));
% st=V(id,1);
% [st,ind]=sort(st,'descend');
% id=id(ind);
% stable=[id st]';
% stable=stable(:);
% stable=stable';
% if length(stable)<6
%     stable=padarray(stable,[0 6-length(stable)],nan,'post');
% else
%     stable=stable(1:6);
% end
% end
% if exist('L_Do','var')
% L_Do=movmean(L_Do,6);
% DDo=L_Do(:,1)+L_Do(:,4)-L_Do(:,2)-L_Do(:,3);
% DDDo=detrend(DDo);
% timp_DDo=1-var(DDDo)/var(DDo);
% DH_Do=detrend(H_Do);
% timp_H_Do=1-var(DH_Do)/var(H_Do);
% 
% temp(soupnumber,:)=[mean(H_Do) timp_H_Do mean(DDo) timp_DDo stable mean(L_Do,1)];
% end
% 
% if exist('NFluxo','var')
% NFluxo(1)=[];
% RFluxo(1)=[];
% ANFluxo(1)=[];
% xs=NFluxo;
% ys=RFluxo;
% zs=ANFluxo;
% xs=movmean(xs,6,'Endpoints','discard');
% ys=movmean(ys,6,'Endpoints','discard');
% zs=movmean(zs,6,'Endpoints','discard');
% dxs=[diff(xs);0];
% dys=[diff(ys);0];
% dzs=[diff(zs);0];
% temp(soupnumber,:)=[mean([xs ys zs],1) mean([dxs dys dzs],1) mean(sum(abs([dxs dys dzs]),2)) estep];
% end
% temp(soupnumber,:)=[mean(E_dmat_V_inputo) std(E_dmat_V_inputo) estep];

% 
% if exist('St1_Ento','var');
% DSt1_Ento=detrend(St1_Ento);
% timp_St1_Ento=1-var(DSt1_Ento)/var(St1_Ento);
% M_St1_Ento=mean(St1_Ento);
% %     St1_Ents=[St1_Ents;M_ST1_Ento timp_St1_Ento];
% 
% DSt2_Ento=detrend(St2_Ento);
% timp_St2_Ento=1-var(DSt2_Ento)/var(St2_Ento);
% M_St2_Ento=mean(St2_Ento);
% %     St2_Ents=[St2_Ents;MSt2_Ento timp_St2_Ento];
% temp(soupnumber,:)=[M_St1_Ento,timp_St1_Ento,M_St2_Ento,timp_St2_Ento,estep];
% %     temp(soupnumber,:)=[mean(coef(:),'omitnan') H_PCAV estep];
% %  temp(soupnumber,:)=[rr,fnrate,mm];
% end
%     

% try
% ind=CD_inputo~=0;
% CD_inputo=reshape(CD_inputo(ind),sum(ind(:))^0.5,[]);
% CCD_inputo=corrcov(CD_inputo);
% CCD_inputo=(CCD_inputo+CCD_inputo')/2;
% [egv eg]=eig(CD_inputo)
% pca();
% eg=diag(eg);
% eg(eg<0)=0;
% eg=eg/sum(eg);
% H_eg=sum(-eg.*log2(eg),'omitnan');
% catch Err
%  CCD_inputo=zeros(18,18);
% H_eg=0;
% end


% egs=[egs padarray(eg,18-length(eg),0,'pre')];
% eso=[eso mean(CCD_inputo(:))];
% 
% a=movmean(H_popo,6);
% b=movmean(H_inputo,6);
% % c=movmean(H_iinputo,6);
% % d=movmean(pt4o,6);
% da=diff(a);
% db=diff(b);
% % dc=diff(c);
% % dd=diff(d);
% % data=[da db dc];
% % dtdata=detrend(data);      
% % timp=1-var(dtdata,1)./var(data,1);
% fa=abs(fft(a));
% fb=abs(fft(b));
% fa(1)=[];
% fb(1)=[];


% hno=hno(31:end);
% popo=popo(31:end);
% input=hno;
% l=length(input);
% window=20;
% hnom=movmean(hno,window);
% mm=max(hno)-min(hno);
% dh=squareform(pdist(hnom'));

% input=hno;
% s=std(input);
% MEAN=mean(input);
% input=(input-repmat(s,1,length(input)))./repmat(MEAN,1,length(input));
% movm=movmean(input',window);
% dh=squareform(pdist(movm));
% dhg=imgaussfilt(dh,0.5);
% rp1=dhg<0.01;
% rp1=rp1(2:end-1,2:end-1);
% rp2=rpmat(dhg,0.001,[1 0 -1; 0 0 0; -1 0 1]);
% dg=conv2(single(diag(diag(rp2))),[0,1,0;1,1,1;0,1,0],'same');
% rp2ndg=rp2 &~ dg;
% 
% mm=max(input)-min(input);
% fnrate=sum(sum(rp1&~rp2))/sum(rp1(:));
% rr=mean(rp2ndg(:));

% [rr,det,entr,l1] = Recu_RQA(dm>0.001,0);

    if update==1;
 
    end
     
str=num2str(temp);
disp(str)
%record e,mamp,mm estep and rule
short=0;
else
    temp(soupnumber,:)=[nan*ones(1,outnumber-1) estep];
     disp('shortlived soup')
     short=1;
end