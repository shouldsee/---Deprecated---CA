
apc=[];
hapc=[];
spec=[];
l=length(input);
l=min(l,100);
input=input(1:l);


if ~exist('window')
window=20;
end
%% calculate APC profile with a moving window
for i=1:l-window
    wd=input(i:i+window);
    a=autocorr(input(i:i+window),window-1);
%     a=hilbert(input(i:i+window));
%     a=fft(wd);
    apc=[apc;a];
    spec=[spec;abs(fft(a))];
    hapc=[hapc;wentropy(a,'shannon')];
end

aapc=[];
% [p,q]=size(apc);
% for i=1:q
%     a=autocorr(apc(:,i),p-1);
%     aapc=[aapc;a'];
% end
% aapc=xcov(apc,20)';
apc(:,1)=[];
% c=cov(apc);
% ic=pinv(c);
% e=gradient(ic);
% e=gradient(e);
d=pdist(apc);
d=squareform(d);


% subplot(2,2,2)
% waterfall(apc)
% plot(p);
% waterfall(aapc)
% imagesc(c)
% surface(ic)
% mat_h=size(d,1);
% d1=zeros(mat_h-1,1);
% d1=[];
% for mh=1:mat_h-1
%     for mh2=1:mh
%         d1=[d1;d(mh,mh2) d(mh+1,mh2+1)];
%     end
% end
% scatter(d1(:,1),d1(:,2))
% md1=mean(d1)
% vd1=var(d1)

%% convolution and retrun maps
FIR1=reshape(str2num('010101010'.').',3,3);
FIR2=reshape(str2num('101000101'.').',3,3);
FIR3=reshape(str2num('111101111'.').',3,3);
FIR=(FIR1*0.5+FIR2)/6;
FIR=reshape(str2num('1001'.').',2,2);
FIR=[-1 -1 ;-1 1];
FIR=[1 0 -1; 0 0 0; -1 0 1];
% FIR=[0 0 0; 0 -4 0; 0 0 0]+FIR1;
% FIR=[-1 -1 -1; -1 8 -1; -1 -1 -1];
input=list{k}';

if strcmp(distmethod(listname),'euclidean')
    s=std(input);
    m=mean(input);

    % input=(input-repmat(s,1,length(input)))./repmat(m,1,length(input));
%     input=movmean(input,window);
end
% k=k+1;


dh=squareform(pdist(input,distmethod(listname)));
dhg=imgaussfilt(dh,0.5);


rp1=rpmat(dhg,0.1,[-1 -1 -1; 0 0 0; 1 1 1]);
rp1s=conv2(dh,[-1 -1 -1; 0 0 0; 1 1 1],'same');
rp1s=rp1s(2:end-1,2:end-1);
rp1=dhg<0.01;
rp1=rp1(2:end-1,2:end-1);
dg=conv2(single(diag(diag(rp1))),[0,1,0;1,1,1;0,1,0],'same');
rp1ndg=rp1 &~ dg;

rp2=rpmat(dhg,0.001,[1 0 -1; 0 0 0; -1 0 1]);
rp2s=conv2(dhg,[1 0 -1; 0 0 0; -1 0 1],'same');
rp2s=rp2s(2:end-1,2:end-1);
dg=conv2(single(diag(diag(rp2))),[0,1,0;1,1,1;0,1,0],'same');
rp2ndg=rp2 &~ dg;

fnrate=sum(sum(rp1&~rp2))/sum(rp1(:));


subplot(2,2,1)
% imagesc(rp1s);
% spy(rp1)
imagesc(dhg)
% caxis([0 0.1])
colorbar
title('raw return map of movmean entropy')

subplot(2,2,2)
imagesc(rp2s)
% spy(rp2);
% caxis([-0.05 0.05])
colorbar
title('return map after convolution')

subplot(2,2,3)
spy(rp1)
title('raw map <= 0.01')
subplot(2,2,4)
spy(rp2)
title('abs(conv map)>=0.001')

% imagesc(dh)
% dhc=cov(dm(:),dh(:));
% invd=diag(diag(dhc))^-0.5;
% dhc=invd*dhc*invd;
% axis([0 0.15 0 0.15])
% xlim([0 1.5])
% ylim([0 1.5])
% surface(gradient(ic,2))
% surface(e)
% t=sprintf('sumcov=%1.2f sumconc=%1.2f',sum(c(:)),sum(e(:)));

% rp=dm>0.001;
% spy(rp)
% [X Y]=find(rp2);
% D=[X Y];
% [ind mdist]=knnsearch(D,D,'k',2);
% mmin=mean(mdist(:,2));
% smin=std(mdist(:,2));
% [rr,det,entr,l1] = Recu_RQA(rp2,0);

% scatter(dm(:),dm1(:))

% histo=hist(dm(:),-2:0.02:2)
% xlim([0 0.2])
% ylim([0 0.5])
% xlabel('lag')
% ylabel('window_start')
% zlabel('autocorrelation of input entropy')
% hist(input)
% plot(hapc);
input=diff(input);