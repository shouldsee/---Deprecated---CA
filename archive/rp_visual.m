%% this script is intended for standalone analysis of entropy series
%% load 'example.mat' to run the scripts below

hno_bak=hno;
popo_bak=popo;

%%
l=length(hno);
length(popo)
%% analyser for rpmat (the currently using one)
figure(2)
clist=[];
rrlist=[];
detlist=[];
entrlist=[];
l1list=[];
dhclist=[];
outlist=[];
list=popolist;
m=length(list);
hlist=[];

distmethod=containers.Map;
distmethod('hnolist')='euclidean';
distmethod('popolist')='euclidean';
distmethod('popiolist')='euclidean';
distmethod('pt2olist')='euclidean';
distmethod('pt2iolist')='euclidean';
distmethod('pt4iolist')='euclidean';
distmethod('pt4olist')='euclidean';
distmethod('cellolist')='hamm';

for k=1:2:m
   
%     print(list{k})
    window=6;
    thres1=0.3;
    thres2=0.001;
    

    figure(2)
    listname='hnolist';
    list=eval(listname);
    analyse
    disp(['input' num2str(k)])
    rr1=sum(rp1(:))/size(rp1,1)^2;
    rr2=sum(rp2(:))/size(rp2,1)^2;
    ruleind=floor((k-1)/3)+1;
    fprintf('mmin=%.5f rr=%.5f,rulenumber=%d',rr1,rr2,ruleind);
    suptitle([listname ' rulenumber=' num2str(ruleind) ' ' rulename{ruleind,:}])

   
    figure(3)
    listname='pt2olist';
    list=eval(listname);
    analyse
    disp(['input' num2str(k)])
    rr1=sum(rp1(:))/size(rp1,1)^2;
    rr2=sum(rp2(:))/size(rp2,1)^2;
    ruleind=floor((k-1)/3)+1;
    fprintf('mmin=%.5f rr=%.5f,rulenumber=%d',rr1,rr2,ruleind);
    suptitle([listname ' rulenumber=' num2str(ruleind) ' ' rulename{ruleind,:}])

    figure(4)
    listname='popolist';
    list=eval(listname);
    input=list{k};
    analyse
    disp(['input' num2str(k)])
    rr1=sum(rp1(:))/size(rp1,1)^2;
    rr2=sum(rp2(:))/size(rp2,1)^2;
    ruleind=floor((k-1)/3)+1;
    fprintf('mmin=%.5f rr=%.5f,rulenumber=%d',rr1,rr2,ruleind);
    suptitle([listname ' rulenumber=' num2str(ruleind) ' ' rulename{ruleind,:}])

    figure(5)
    listname='pt4olist';
    list=eval(listname);
    analyse
    disp(['input' num2str(k)])
    rr1=sum(rp1(:))/size(rp1,1)^2;
    rr2=sum(rp2(:))/size(rp2,1)^2;
    ruleind=floor((k-1)/3)+1;
    fprintf('mmin=%.5f rr=%.5f,rulenumber=%d',rr1,rr2,ruleind);
    suptitle([listname ' rulenumber=' num2str(ruleind) ' ' rulename{ruleind,:}])

% input=cumsum(input);
% dhc(2,1)
% pause
% rp2=crp(movm,'nrmnorm','silent','true');
% fn=fnn(input,'silent',true);
% % pause
% mat3=conv2(single(dh),[1,0,1;0,0,0;1,0,1],'same');
% mat3=mat3(2:end-1,2:end-1);
% mat3=mat3/sum(mat3(:));
% mat4=conv2(single(dh),[0,1,0;1,0,1;0,1,0],'same');
% mat4=mat4(2:end-1,2:end-1);
% mat4=mat4/sum(mat4(:));
% 
% % dg=conv2(single(diag(diag(rp2))),[0,1,0;1,1,1;0,1,0],'same');
% rp2ndg=rp2 &~ dg;
% h=hist(mat3(rp2ndg),-0.2:0.01:0.2);
% nh=h/sum(h);
outlist=[outlist; rr1 rr2 fnrate];
% hlist=[hlist;nh];

% clist=[clist;c(2,1)];
% end


% figure(7)

% 
% clf
% subplot(2,2,1)
% % imagesc(d)
% imagesc(mat4.*rp2ndg)
% colorbar
% % t=sprintf('A complex rule for 500 steps with window of 20 \n');
% % title(t);
% subplot(2,2,2)
% imagesc(dh);
% colorbar
% subplot(2,2,3)
% spy(rp2)
% subplot(2,2,4)
% imagesc(abs(mat3).*(rp2ndg))
% % stem(h)
% colorbar
pause
end
s=std(outlist,1);
s=repmat(s,length(outlist),1);
outlist=outlist./s;
% waterfall([rrlist/sum(rrlist);detlist/sum(detlist);entrlist/sum(entrlist);l1list/sum(l1list) dhclist/sum(dhclist)]);
% subplot(2,2,4)
% imagesc(hlist)
subplot(1,1,1)
waterfall(outlist')
zlim([-2 8])

% plot(outlist)

%%

clf
list=popolist;
m=length(list);
timpms=[];
for rind=1:floor((m)/soupmax)+1
    figure(8)
    fprintf('mmin=%.5f rr=%.5f,rulenumber=%d',rr1,rr2,rind);
    suptitle([listname ' rulenumber='  num2str(rind) ' ' rulename{rind,:}])
    cellos=[];
    cs=zeros(4,4);
    cs=[];
    bs=[];
    timps=[];
    for k=(soupmax*(rind-1)+1):soupmax*rind
%         for k=soupmax*(rind-1)+1
    S_pop=S_poplist{k}';
    if sum(S_pop(:))~=0
        x=movmean(H_poplist{k},6);
        y=movmean(H_input{k},6);
        z=movmean(H_iinput{k},6);
%         a=movmean({k},6);
        t=1:length(x);
        
        dx=diff(x);
        dy=diff(y);
        dz=diff(z);
        da=diff(a);
        dt=diff(t);
        subplot(2,3,1)
%        scatter3(x,y,z,'Linewidth',1)
%         quiver3(x,y,z,[dx 0],[dy 0],[dz 0],'MaxHeadSize',5)
%         text(x(1),y(1),z(1),num2str(k))
        quiver(x,y,[dx 0],[dy 0],'MaxHeadSize',5)
        text(x(1),y(1),num2str(k))
%         line(x,y,z,'Linewidth',1)
        hold on
        subplot(2,3,2)
%         quiver3(x,y,z,[dx 0],[dy 0],[dz 0],'MaxHeadSize',5)
%         text(x(1),y(1),z(1),num2str(k))
        quiver(x,y,[dx 0],[dy 0],'MaxHeadSize',5)
        text(x(1),y(1),num2str(k))
%       
%         scatter3(dx,[dx(2:end) 0],(3*k)*ones(1,length(dx)));
%         scatter3(dy,[dy(2:end) 0],(3*k+1)*ones(1,length(dx)));
%         scatter3(dz,[dz(2:end) 0],(3*k+2)*ones(1,length(dx)));
        
        hold on
                
        cellos=[cellos;S_pop(1:40,:)];
        dxt=dx/sum(dx);
        dxh=sum(-dxt.*log2(dxt),'omitnan');
        
%         hist(dxt);
         c=cov([dx;dy;dz;da;t(1:end-1)]');
        invc=diag(diag(c))^-0.5;
        c=invc*c*invc;
        
%         cs=(cs+c);
        cs=[cs;c];
        b=(xcorr([x;y;z]'));
        bs=[bs;b];
        fprintf('\n %d ',k)
        fprintf(['dxh=' num2str(mean(dx))])

        data=[dx;dy;dz;da]';
        dtdata=detrend(data);
        
        timp=1-var(dtdata,1)./var(data,1);
        timps=[timps;timp];
         end
        
    end
    timpms=[timpms;mean(timp,1)];
    subplot(2,3,1)
    hold off
    axis([0 2 0 3 0 4])
    axis([0 1 0 3])
    subplot(2,3,2)
%     axis([-0.1 0.1 -0.1 0.1])
    hold off
    subplot(2,3,3)
    imagesc(cs)
    subplot(2,1,2)
%     waterfall(timps')
    zlim([0 0.5])
    sumcov=mean(cs(:),'omitnan');
    title(sprintf('sum-cov=%f',sumcov))
    caxis([0 1])
    colorbar
    
%     figure(10)
%     [cr,lgs] = xcorr([dx;dy;dz]','coeff');

%     for row = 1:3
%         for col = 1:3
%             nm = 3*(row-1)+col;
%             subplot(3,3,nm)
%             stem(lgs,cr(:,nm),'.')
%             title(sprintf('c_{%d%d}',row,col))
%             ylim([-1 1])
%         end
%     end
%     subplot(2,3,4);
%     autocorr(dx)
%     subplot(2,3,5)
%     autocorr(dy)
%     subplot(2,3,6)
%     autocorr(dz)
%      spy(cellos);
%     imagesc(cellos);
%     d=squareform(pdist([x;y;z]'));
%     subplot(2,3,2)
%     imagesc(d)
%     imagesc(bs)

    
    figure (9)
clf
    cellos=[];
    cs=zeros(4,4);
    cs=[];
    suptitle([listname ' rulenumber='  num2str(rind) ' ' rulename{rind,:}])
    
    for k=(soupmax*(rind-1)+1):soupmax*rind
%         for k=soupmax*(rind-1)+1
 
        S_pop=cellolist{k}';
        if sum(S_pop(:))~=0
            cellos=[cellos;S_pop];
        end
         
        c=cov([x;y;z;a]');
        invc=diag(diag(c))^-0.5;
        c=invc*c*invc;
        
%         cs=(cs+c);
        cs=[cs;c];
        disp(k)
    end
  imagesc(cellos);  
    pause
end

%%
figure(11)
waterfall(timpms);
%%
subplot(1,3,1)
d=[input];
f=[abs(fft(input))];
r=[input];
l=length(input);
plot(input)
hold on
for i=1:5
    dt=diff(input,i);
    rc=cumsum(dt,i);
    d=[d;dt zeros(1,l-length(dt))];
    r=[r;rc zeros(1,l-length(rc))];  
    f=[f;abs(fft(dt)) zeros(1,l-length(dt))];
    plot(rc)
end
hold off

figure(8)
subplot(1,1,1)
waterfall(d);
% subplot(1,3,3)
% waterfall(r);
%%
a=autocorr((hno),min(160,length(hno)-1));
 b=autocorr(popo,min(160,length(popo)-1));
Fs = 1;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = length(a);             % Length of signal
t = (0:L-1)*T;        % Time vector
Y = fft(a);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure(9)
subplot(3,1,1)
autocorr(hno,min(160,length(hno)-1));
subplot(3,1,2)
plot(f,(P1))
title(  sprintf('entropy of spectrum=%f \n, max amp=%f, sum amp=%f',we,maxamp,sum(P1)))
xlabel('f (Hz)')
ylabel('|P1(f)|')
subplot(3,1,3)
Y = fft(b);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,(P1))
%autocorr(a,min(160,length(a)-10));
%         plot(sort(P1,'descend'))
title(sprintf('ruleB=[%s] mm=%f',num2str(ruleB),mm));