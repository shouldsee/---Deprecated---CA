global cells cells_old n x y rulecurr async
n=200;cooldown=100;
x=2:n+1;
y=2:n+1;
cells=zeros(n+2,n+2);
async=0;
pw=3;
rsvct=shiftdim(2.^(0:pw-1),-1);
async=0;
siz=size(cells);
[x1,x2]=ndgrid(x,y);
xyid=sub2ind(siz,x1,x2);
n3=floor(n/3);

ix=n3:2*n3;
cellid=1:numel(cells);
border=cellid(~ismember(cellid,xyid));
FIRM=ones(3,3)/9;
% ix=x(~ismember(x,ix));
[ix1,ix2]=ndgrid(ix,ix);
ixid=sub2ind(siz,ix1,ix2);
oxid=xyid(~ismember(xyid,ixid));
lag=pw^2+pw-1;

lag2=20;ix=n3-10:2*n3+10;
[ix1,ix2]=ndgrid(ix,ix);
ixid=sub2ind(siz,ix1,ix2);

names={'nerd','weirdo','visual_chaos'};
randrule=0;
for i=1:numel(names);
	name=names{i};
		
    if strcmp(name,'rand');
		randrule=1;
		else
			randrule=0;
		end
		loadrule
		rind=1;
		lgc;
end
%%

name='rand';
loadrule;
rulename=repmat(rulename(1),2E6,1);
% rulename=rulename';
randrule=1;
rind=1;
lgc;
%%
name='defnt';
% loadrule
load(name,'rulename')
rulename=repmat(rulename(1),2^17,1);
rulecurr=zeros(1,18);
randrule=2;
rind=1;
lgc
%%
name='defnt2';
% loadrule
load(name,'rulename')
rulename=repmat(rulename(1),2^17,1);
rulecurr=zeros(1,18);
randrule=2;
rind=2^17+1;
lgc


%%
names={'nerd','weirdo','visual_chaos','rand'};
for i=1:numel(names);
    name=names{i};
    loadrule
    rs.(name)=rulename;
    
end
%%
%%
rec=ak.defnt;
% recs=rec(:,1:end-3);
pws=rec(:,end-2:end);
ind=find(pws(:,1)<3.75 & pws(:,2)>0.15 & pws(:,2)<0.5 & pws(:,3)<1.9);
% rec=rec(ind,:);
rs.defntf=rs.defnt(ind);
ak.defntf=rec(ind,:);

pws=[];
rules=[];
recs=[];
group=[];


names={'nerd','weirdo','visual_chaos','rand'};
names={'defntf','weirdo','visual_chaos','nerd'};
for i=1:numel(names);
	
		name=names{i};
% 		loadrule
		rec=ak.(name);
		recs=[recs;rec(:,1:end-3)];
		pws=[pws;rec(:,end-2:end)];
        rulename=rs.(name);
        rulename=rulename(1:size(rec,1));
		rules=[rules;rulename];
		group=[group;i*ones(size(rec,1),1)];

end

figure(3)
clf
fig=gcf;
% colormap colorcube
scatter3(pws(:,1),pws(:,2),pws(:,end),10,group);
fig=gcf;
dcm_obj=datacursormode(fig);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,rules});
