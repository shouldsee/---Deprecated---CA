%%
inds=(1:length(rec0))';
i=1;
rec=rec0;
while i<=length(rec);
d1=pdist2(rec(i,:),rec);
siminds=(d1>0.5);
rec=rec(siminds,:);
inds=inds(siminds,:);
% inds=inds(~ismember(inds,siminds));
% inds(siminds)=[];
i=i+1;
fprintf('i=%d, matl=%d \n',i,length(rec));
end
% rind=rind+1;
%%
rss=[];
r0=rulecurr;
for i=1:18
    r=r0;
    r(i)=1-r(i);
    r=(base2dec(num2str(r')',2));
    rss(i,:)=r;
end
%%
for i=1:length(id)
rulecurr=str2num(dec2base(id(i)-1,2,18)')';
    rss(i,:)=rulecurr;
end
%%
[~,wh]=ismember([0 0 1 0 0 0 0 0 1 0 0 0 0 1 1 1 1 1],rss,'rows')
%%
figure(3)
clf
dt=squareform(pdist(rss,'hamm'));
[V D]=eigs(dt,2);
[V,iid]=sort(V(:,2));
imagesc(dt(iid,iid))
%%
lst=find(dt(rind0,:)<=1/18)
axis equal
%%
clf
plot(dt(,:))
%%
d=squareform(pdist(rec(1:2822,:),@fdist));
%%
Z=linkage(rec,'average');
c=cophenet(Z,d);
%%
[V D]=eigs(d,3,'SA');
[V,ind]=sort(V(:,2));
imagesc(d(ind,ind));
%%
figure(5)
id=inds(ind);
imagesc(rec(ind,:));