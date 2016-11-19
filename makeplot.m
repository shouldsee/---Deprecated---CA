if ~exist('figh')
    figh=5;
end
% figure(figh)
hold on
% histogram(pdata,0:0.025:4)
xlim([0 4])

%%
h=figure(4);
% set(gcf, 'Visible', 'off')

% Hxy=H_inputo(2:end-1,2);
% Hyz=H_inputo(3:end,2);
% Hxyz=H_inputo(3:end,3);
% Hy=H_inputo(2:end-1,1);
% Hz=H_inputo(3:end,1);
% Iinfo=Hxy+Hyz-Hxyz-Hy;
input=S_input;
CS_inputX=circenum(input,[1 0]);
CS_inputY=circenum(input,[0 1]);
CS_inputX=reshape(CS_inputX,size(CS_inputX,1),[]);
CS_inputY=reshape(CS_inputY,size(CS_inputY,1),[]);
dCS_input=pdist2(CS_inputX,CS_inputY,@calc_MI);
dCS_input=circshift(full(dCS_input),[15 15]);
dCS_input=dCS_input/dCS_input(16,16);
ldCS_input=real((dCS_input));


% 
% subplot(3,2,1);
% input=H_inputo;
% plot(input*var(input))
% hold on
% input=movmean(input,6,'Endpoints','discard');
% plot(input*var(input))
% hold on
% apc=autocorr(input,min(142,length(input)-1))*var(input);
% var(input)
% apc(mod(1:length(apc),6)~=1)=0;
% plot(apc*4);
% hold off
% ylim([-1 1]*0.05)


% subplot(3,2,2);
% input=MI_vico;
% apc=autocorr(input,min(50,length(input)));
% plot(apc);
% hold on
% apc=movmean(apc,6,'Endpoints','discard');
% plot(apc);
% hold on
% plot(input);
% hold off
% % ylim([-1 1])

% imagesc(log2(D_defecto));
% caxis([0 5])
% plot(MI_inputo)
% plot(ldCS_input(15,:));
%  xlim([0 30]);
% imagesc(ldCS_input);
% ylim([-5 1]);
% plot(KL_divo)
% plot(H_inputo(:,1),HS_vico);
% plot(HS_vico);
% xlim([0 140]);
% ylim([0 0.5]);
% hold on
% plot(H_input_oldo);

% plot(MI_vico/H_input_oldo);
% MI_vico=2*H_input_oldo-HS_vic2o;
% xs=H_popo;
% ys=H_inputo;
% zs=H_dinputo;
% so=[so;MIso rulename(rind,:)];
% ruleo=
% xs=HS_vic2o;
% xs=0.5*(H_dinputo(1:end-1)+H_dinputo(2:end));
% ys=MI_dinputo(2:end);
% ys=HSS_vic2o;
% xs=H_dinputo;
% ys=H_dinput_inputo;
% xs=H_input_oldo;
% ys=MI_vico;
% xs=diff(xs);
% ys=diff(ys);
% scatter3(xs,ys,1:length(xs),3,'x');
% scatter3(xs,ys,zs,3,'x');
% plot(MIso)
% scatter(xs,ys,3,'x');

% c=corrcov(cov(xs,ys));
% title(['correlation' num2str(regress(ys,xs))]);
% xlim([0 10])
% ylim([0 10])
% zlim([0 10])
% subplot(3,2,2);
% FS_flowo=bsxfun(@(a,b) a.*(a>b),LS_flowo,prctile(LS_flowo,95,2));
% FS_flowo(FS_flowo==0)=nan;
% MFS_flowo=mean(FS_flowo,2,'omitnan');
% M_flowo=mean((LS_flowo-1).^2,2).^0.5;
% P_flowo=bsxfun(@rdivide,D_flowo,sum(D_flowo,2));
% SK_flowo= skewness(LS_flowo,[],2);
% HD_flowo=log2(D_flowo);
% HD_flowo(isnan(HD_flowo))=0;
% 
% KL_flowo=kurtosis(log2(LS_flowo),[],2);
% HP_flowo=-P_flowo.*log2(P_flowo);
% HP_flowo(isnan(HP_flowo))=0;
% % H_flowo=sum(HP_flowo,2);
% % plot(H_inputo,HLEo);
% xlim([0,4])
% ylim([0 14])
% ylim([0 1] * 15)
% xs=H_inputo;
% ys=FH_inputo;
% ys=N_defecto;

% xs=2.^(-H_inputo);

% ys=MI_popo;
% xs=H_input_oldo;
% ys=MI_vico;
 % xs=diff(xs);
% ys=diff(ys);
% lxs=makezero(-log2(xs));
% lys=makezero(-log2(ys));

% scatter(lxs,lys,13,'x');
% [coeff,~,latent,~,explained]=pca([lxs,lys]);
% wcoeff=bsxfun(@times,latent,coeff');

% plot(xs,ys);
% scatter3(lxs,lys,1:length(lxs),3,'x');

% plot3(lxs,lys,1:length(lxs));
% scatter(xs,ys,3,'x');
% xlim([-1 1]*0.2)
% ylim([-1 1]*0.2)
% xlim([0,4])
% ylim([0,4])
% plot(ys./xs)

% 
% dLS_popo=diffcomb(LS_popo);
% dLS_popo=dLS_popo(1:2:end,:);
% % dmat=squareform(pdist(LS_popo,@calc_MI));
% dmat=real(log2(dCS_input));
% imagesc(dmat)
% caxis([-4 0])

subplot(3,2,3);
% plot((N_defecto./(2.^(-H_inputo))));
xlim([0 150]);
% xlim([0 30])
% ylim([-5 0])
ylim([0 20])
% ylim([0 250]);
% maxH=@(a,b) max([calc_Hx(a,b),calc_Hy(a,b)],[],2);
% ddLS_popo=diffcomb(dLS_popo);
% dmat=squareform(pdist(LS_linputo,@calc_MI));
% dmat=real(log2(dmat));
% imagesc(HSD_vic)

subplot(3,2,4)
% plot(movmean(abs(diff(N_defecto./(2.^(-H_inputo)))),6,'Endpoints','discard'))
% dmat=squareform(pdist(LS_inputo,@calc_MI));
% dmat=real(log2(dmat));
% imagesc(dmat)
% ami3
% scatter3(H_inputo,Hxy2o,HH_Tinputo,2,'x');
xlim([0 150])
% ylim([0 8])
% zlim([0 8])

m=size(reco,1);

% pdata=reshape(reco(2:1+m-mod(m,warmup),:,:),[],warmup,warmup,3);
vreco=squeeze(reshape(reco,[],1,3));
gp=findgroups(vreco(:,2),vreco(:,1));
mvreco=splitapply(@mean,vreco(:,3),gp);
mvreco=reshape(mvreco,[],warmup)';
% mvreco=flip(mvreco,1);
% mvreco=flip(mvreco,2);
% 
subplot(2,2,1)
% mpdata=squeeze(mean(pdata,1));
imagesc(mvreco');
view(0,90);
set(gca,'YTick',1:length(po));
a=po;
yticks=cellfun(@num2str,num2cell(a),'UniformOutput',false);
set(gca,'YTickLabels',yticks);
ylabel('Initial Density')
xlabel('Nth step from randomness')
title('mean(H input)')

colorbar

subplot(2,2,2)
svreco=splitapply(@std,vreco(:,3),gp);
svreco=reshape(svreco,[],warmup)';
% spdata=squeeze(std(pdata,1));
imagesc(svreco');

a=po;

xlabel('Nth step from randomness')
title('std(H input)')
colorbar

subplot(2,2,3)
msvreco=svreco./mvreco;

imagesc(msvreco');

set(gca,'YTick',1:warmup);

a=po;
yticks=cellfun(@num2str,num2cell(a),'UniformOutput',false);
set(gca,'YTickLabels',yticks);

xlabel('Nth step from randomness')
title('std(H input)')
colorbar


suptitle(['MISC' ' rulenumber='  num2str(rind) ' ' rulename{rind}])

% subplot(3,2,[3 4 5 6])
subplot(2,2,[4])
% figure(5)
% clf
xlim([0 1])
% a=splitapply(@mean,reco(:,2),reco(:,1)+1)
% hold on
% pdata=squeeze(reshape(reco,1,[],3));
% pdata(:,3)=squeeze(std(reco(:,:,3),1));
% pdata(:,[1 2])=squeeze(mean(reco(:,:,[1 2],1)));
% pdata
% % ppdata=reshape(pdata,[],3);
ppdata=reshape(reco,[],3);
scatter3(ppdata(:,1),ppdata(:,2),ppdata(:,3),3,'x')
hold on
surface(po,1:warmup,mvreco);
hold off

% zlim([0 2])
ylim([1 warmup])
view(60,20)
% histogram2(reco(:,1),reco(:,2),...
%     [max(reco(:,1)),100],...
%     'Normalization','pdf')

% view(-90,90)

% xlim([0 max(reco(:,1))])
% ylim([0 6])
% zlim([0 1])

% imagesc(LSI_inputo);
% subplot(3,2,[5 6])
ylabel('Nth step from randomness')
xlabel('initial density')
zlabel('H input')
% title('1D-ised space-time pattern')
% subplot(3,2,[3 4])
% imagesc(LCS_inputo);
% caxis([0 1])
% title('1D-ised space-time pattern')
drawnow 
if store==1;
    makeprint
end
