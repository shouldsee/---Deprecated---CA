
%%
if ~exist('randrule','var')
    randrule=0;
end
if ~exist('randcell','var')
    randcell=1;
end
if ~exist('record','var')
    record=1;
end
if ~exist('gifname','var')
    gifname='gif';
end

global cells cells_old n x y rulecurr async
rind=1;
k=rind;
% randrule=1;
% randrule=0;
name='nerd';
% loadrule
if randrule
getrule;
end
n=300;
x=2:n+1;
y=2:n+1;
if randcell
cells=zeros(n+2,n+2);
end
async=0;
siz=size(cells);
[x1,x2]=ndgrid(x,y);
xyid=sub2ind(siz,x1,x2);
p0=0.5;

%%
if randcell
cells=rand(n+2,n+2)*1;
% cells=zeros(n+2,n+2);
cells(10:15,10:15)=rand(6,6)/2;
end
%
cells=torus(cells);
rulecurr=@(a) 1-(mod(100*a,100)/100);
% rulecurr=@(a) (mod(a,4))/4;
rulecurr=@(a) tanh(a/10).*sin(a)

figure(1)
fi=imagesc(cells(xyid)',[-0.5 1]*1);
figure(2)
h=histogram(cells(xyid));
figure(3)
h2=histogram2(cells(xyid),cells(xyid),100);
edge=-0.5:0.025:1;
h2=histogram2(cells(xyid),cells(xyid),edge,edge);
view(30,30);
zlim([10,2000]);
set(gca,'ZScale','log')
ex=15;
px=linspace(-ex,ex,n+2);
py=linspace(-ex,ex,n+2);
px=px(x1);
py=py(x2);
if ~exist('stepmax','var')
    stepmax=150;
    intl=2;
end
    stepnum=1;

div=1E80;
f1='Sinput(xyid)=mod(Sinput(xyid),px)./py;';
f2='Sinput(xyid)=cos(Sinput(xyid)./px)./py;';
ef1='eval(f1);';
ef2='eval(f2);';

% updateFcn=['Sinput(xyid)=eval(f1)./eval(f2);'];
% updateFcn=['Sinput(xyid)']
% 'cellsc=cells(xyid)'

if randrule
[ux,uy]=meshgrid(-1:9,-1:2);

% uz0=uz;
    uz=((uy==0 & ismember(ux,B) )| uy==1 & ismember(ux,S));
% uz=rand(size(uz))>0.5;
uz=single(torus(uz));
upf=@(a,b) interp2(ux,uy,uz,a,b);
end



upf=@(a,b) mod(a,8)+mod(b,1);
fc='Sinput(xyid)=arrayfun(upf,Sinput(xyid),cells(xyid));';
fc='Sinput(xyid)=mod(mod(Sinput(xyid),8)./px+mod(cells(xyid),1)./py,1);';
% fc='Sinput(xyid)=(mod(Sinput(xyid),8)./px+mod(cells(xyid),1)./py)./(8./px+1./py);';
fc='Sinput(xyid)=(mod(Sinput(xyid),8)./px+mod(cells(xyid),1)./py);';
updateFcn=[fc];
% f2='Sinput(xyid);'

% px=px(128,1);
% py=py(1,96);
% n=60; ex=500;
% rect=[152,313;
%       214,350];
% rect=[365,80;
%       465,135];
%   rect=[514,238;
%       541,252];
  rect=[207,301;
      256,312];
rect=[281,288;
      320,320];
rect=[306,302;
      350,315];
rect=[1,1;
      298,298];
  rect=[132,325;
      150,326];
rect=[271,1;
   500,231];
% rect=[80,80;160,160];
% % rect=[350,100;500,250];

% rect=[40,80;
%       160,81];
% rect=[1,20;200,150];
rect=[300,140;329,155];
rect=[163,60;164,60];
px=repmat(linspace(px(rect(1,1),1),px(rect(2,1),1),n)',1,n);
py=repmat(linspace(py(1,rect(1,2)),py(1,rect(2,2)),n),n,1);
% rect=[193,209;194,210];
% rect=[200,96;200,96];
% rect=[51,64;126,121];
% px=repmat(linspace(px(rect(1,1),1),px(rect(2,1),1),n)',1,n);
% py=repmat(linspace(py(1,rect(1,2)),py(1,rect(2,2)),n),n,1);

% 
% rect=[70,93;
%       77,96];
% % rect=[73,93;73,93];
% px=repmat(linspace(px(rect(1,1),1),px(rect(2,1),1),n)',1,n);
% py=repmat(linspace(py(1,rect(1,2)),py(1,rect(2,2)),n),n,1);
% 
% 
% px=repmat(linspace(px(25,1),px(35,1),n)',1,n);
% py=repmat(linspace(py(1,135),py(1,150),n),n,1);

xlabel('px');
ylabel('py');
ax=fi.Parent;
xlabel(ax,'px');
ylabel(ax,'py');
set(ax,'XTickLabels',cellstr(num2str(px(ax.XTick,1),3)));
set(ax,'YTickLabels',cellstr(num2str(py(1,ax.YTick)',3)));
set(ax,'XTickLabels',cellstr(num2str(px(ax.XTick,1),3)));
set(ax,'YTickLabels',cellstr(num2str(py(1,ax.YTick)',3)));
% cells=gpuArray(cells);
mvs=zeros(1,stepnum);
mv=mean(cells(:));
ck0=eye(n,n);

for stepnum=1:stepmax
    %%
    pr=rand(n,n);
cells=floor(div*cells)/div;
cold=gather(cells(xyid));

S_input=conv2(cells,[1 1 1; 1 0 1; 1 1 1],'same');
Sinput=(floor(div*S_input)/div);
% cells=rulecurr(Sinput);
% Sinput(xyid)=eval(f1);
% Sinput(xyid)=eval(f2);
eval(updateFcn);
% Sinput(xyid)=mod(Sinput(xyid),px)./py;
cells(xyid)=Sinput(xyid);
cells=torus(cells);
cellsT=gather(cells(xyid));
%stepnum=stepnum+1;
if mod(stepnum,intl)==0
set(fi,'CData',cells(xyid)')
set(h,'Data',cells(xyid))
set(h2,'Data',[cold(:),cellsT(:)])
mv=mean(cells(:));
MAX=max(cells(:));
MIN=min(cells(:));
tl=sprintf([
    'stepcount=%d \n',...
    'div=%.2d,\n',...
    'mean=%.2d ' ,...
    'min=%.2d max=%.2d \n',...
    '%s \n'],stepnum,div,mv,MIN,MAX,updateFcn);
% set(gca,'Title','1')
title(gca,tl)
drawnow
if record
im=frame2im(getframe(fi.Parent.Parent));
[imind,cm]=rgb2ind(im,256);

if stepnum==intl;
    imwrite(imind,cm,[gifname '.gif'],'gif','Loopcount',inf,'DelayTime',0.05);
else
    imwrite(imind,cm,[gifname '.gif'],'gif','WriteMode','append','DelayTime',0.05);
end
end
pause(0.05)
% stepnum=stepnum-intl+1;
end
mvs(stepnum)=gather(mv);
end

%%
figure(6)
histogram(cells(xyid),linspace(0,0.4,100));
figure(7);
plot(linspace(0,1,length(mvs)),abs(fft(mvs)))
%%
max
(px(:)); min(px(:))