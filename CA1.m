
init_CA
getrule;

pw=2;
lag=2^pw+pw-1;
cells=single(rand(n+2,n+2)<0.5);
figure(1)
subplot(1,1,1)
fi=imagesc(cells(xyid)',[0 1]*1);
figure(2)
subplot(2,1,1)
h2=histogram2(cells(xyid),cells(xyid),100);
set(gca,'ZScale','log')
subplot(2,1,2)
h1=histogram(cells(xyid),100);
set(gca,'YScale','log')
% subplot(2,2,2)
% h1=histogram(1:2^20,4000);
% set(gca,'YScale','log')
% subplot(2,2,3)
stepmax=1000;
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

for stepnum=1:stepmax
   Sinput=conv2(cells,FIR.S_input,'same');
   cells=rulecurr(Sinput+1);
   cells=torus(cells);
   ci=reshape(cells(xyid),1,[]);
   cold=[ci;
       cold(1:lag-1,:)];
   Eg1T=edge(H_cpdT,'sobel');
   Eg1T=torus(Eg1T);
   Eg1Tstk=cat(3,Eg1T,Eg1Tstk(:,:,1:lag-1)); 
   if stepnum>lag
    
    cpd=conv2(cold,2.^(0:pw-1)','valid');
    ct_cpd=sum(hist(cpd,0:2^pw-1)~=0,1);
    H_cpd_old=H_cpd;
    H_cpd=reshape(log(ct_cpd)/log(2^pw),n,n);
    H_cpdT(xyid)=H_cpd;
    H_cpdT=torus(H_cpdT);
%     H_cpdTstk=cat(3,H_cpdT,H_cpdTstk(:,:,1:lag-1));
%     H_cpdTsmt=mean(H_cpdTstk,3);
    
    cH_cpdT=conv2(H_cpdT,[1 1 1; 1 1 1;1 1 1]/9,'same');
    cH_cpd_old=cH_cpd;
    cH_cpd=cH_cpdT(xyid);
    Eg1Tsmt=mean(Eg1Tstk,3)>0.5;
    fEg1T=conv2(single(Eg1T),[0 1 0;1 1 1;0 1 0],'same')~=0;
    Eg1=fEg1T(xyid);
    
%     imagesc(H_cpd.*Eg1)
%        rsd_cold=bsxfun(@times,cold,expos);
%        crsd_cold=sum(rsd_cold(1:lag-1,:),1);
%    HDD=hist(crsd_cold,(1:2^20)-1);
%    HD=entropise(HDD);
%    H=sum(HD);
%    nMI=calc_MI(crsd_cold,reshape(mod(Sinput(xyid)-1,9),1,[]))
%    rec(stepnum-lag)=nMI;
   end
   if mod(stepnum,intl)==0
%       set(fi,'CData',cells(xyid));
%       set(fi,'CData',H_cpd.*(abs(cH_cpd))); 
      set(fi,'CData',H_cpd); 
       set(h2,'Data',[cH_cpd_old(:),cH_cpd(:)]);
    set(h1,'Data',H_cpd)
   end
%    set()
  drawnow
     
end
subplot(2,2,3)   
% rec=movmean(rec,2,'Endpoints','Discard');
% plot(rec,[nan diff(rec)]);
rulename{rind}
rind=rind+1

% ylim([-0.05 0])