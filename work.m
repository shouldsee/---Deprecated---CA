load visual_chaos
scanthrough
col=1;
EDVI=E_dmat_V_inputos;
gp1=col*ones(1,length(E_dmat_V_inputos))';
load visual_complex
scanthrough
EDVI=[EDVI;E_dmat_V_inputos];
col=col+1;
gp1=[gp1;col*ones(1,length(E_dmat_V_inputos))'];

%%
S=S_poplist;
R=rulename;
load('visual_ambig','S_poplist','rulename')

S_poplist=[S_poplist;S];
rulename=[rulename;R];
%%
load visual_ambig
scanthrough
col=col+1;
EDVI=[EDVI;E_dmat_V_inputos];
gp1=[gp1;col*ones(1,length(E_dmat_V_inputos))'];
%%
% gp1(end-length(E_dmat_V_inputos)+1:end)=[];
l=length(E_dmat_V_inputos);

%%
ind=length(EDVI)-l+1:length(EDVI);
%%
% ind=1:20000;
pdata=EDVI(ind,:);
cols=gp1(ind);
%%
figure(7)

subplot(2,2,1)
xs=pdata(:,1);
ys=pdata(:,2);
zs=pdata(:,3);
cs=cols;
scatter3(xs,ys,zs,3,cs*10,'x');

subplot(2,2,2)
xs=pdata(:,4);
ys=pdata(:,5);
zs=pdata(:,6);
cs=cols;
% scatter3(xs,ys,zs,3,cs*10,'x');
scatter(ys,zs)