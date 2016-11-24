
init_CA
getrule;

pw=4;
dm=2;
% lag=max(ceil(2^(pw^dm)/n^2),2);
lag=2;
cells=single(rand(n+2,n+2)<0.5);
figure(1)
subplot(1,1,1)
fi=imagesc(cells(xyid)',[0 1]*1);
figure(2)
% subplot(2,1,1)
% h2=histogram2(cells(xyid),cells(xyid),100);
% set(gca,'ZScale','log')
% subplot(2,1,2)
% h1=histogram(cells(xyid),100);
% set(gca,'YScale','log')
% subplot(2,2,2)
% h1=histogram(1:2^20,4000);
% set(gca,'YScale','log')
% subplot(2,2,3)
cold=zeros(lag,n^2);
expos=2.^((1:lag)-1)';
rec=[];
H_cpd=cells(xyid);
H_cpd_old=H_cpd;
H_cpdT=cells;
cH_cpd=cells(xyid);
cH_cpd_old=cells(xyid);
Eg1Tstk=repmat(cells,1,1,lag);
H_cpdTstk=repmat(cells,1,1,lag);
coldT=repmat(cells,1,1,lag);
RsVct=2.^(reshape(0:pw^dm-1,repmat(pw,1,dm)));
vnum=4;

rulemax=10000;
reco=zeros(rulemax,stepmax-lag);
expos=2.^((1:lag)-1)';

%%
for rulenum=rulenum0:rulemax;
cells=randi([0 1],[n+2 n+2]);
rec=zeros(1,stepmax-lag);
rind=rulenum;
getrule
for stepnum=1:stepmax
   Sinput=conv2(cells,FIR.S_input,'same');
   cells=rulecurr(Sinput+1);
   cells=torus(cells);
   coldT=cat(3,cells,coldT(:,:,1:lag-1));
   order=spatial_corr(cells);
%    rsd_coldT=convn(coldT,RsVct,'same');
%    rsd_cold=rsd_coldT(x,y,:);
%    ct=sum(hist(rsd_cold(:),0:2^(pw^dm)-1)~=0);
%    H_cpd=log(ct/2^(pw^dm));
%    order=-H_cpd;
   if stepnum>lag
   rec(stepnum-lag)=order;
   if stepnum>vnum+lag;
   if sum(diff(rec(stepnum-lag-vnum+1:stepnum-lag))<0.005)==vnum;
       break
   end
   end
   end
   if mod(stepnum,intl)==0
      set(fi,'CData',cells(xyid));
%       set(fi,'CData',H_cpd.*(abs(cH_cpd))); 
%       set(fi,'CData',H_cpd); 
%        set(h2,'Data',[cH_cpd_old(:),cH_cpd(:)]);
%     set(h1,'Data',H_cpd)
   end
%    set()
  drawnow
     
end
reco(rulenum,:)=rec;

end
%%
figure(2)
% hold on
rec=movmean(rec,2,'Endpoints','Discard');
plot(rec,[nan diff(rec)]);
hold off
rulename{rind}
rind=rind+1
% order
xlim([0 5])
ylim([0 0.2])
%%
subreco=reco(1:rulenum,:);
[wcoeff,score,latent,tsquared,explained]=pca(subreco);
%%
c3 = wcoeff(:,1:3);
coefforth = inv(diag(std(subreco')))*wcoeff;
I = c3'*c3;
cscores = zscore(subreco')*coefforth;
%%
figure(3)
dcm_obj = datacursormode(3);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,rulename})
plot3(score(:,1),score(:,2),score(:,3),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
% ylim([-0.05 0])