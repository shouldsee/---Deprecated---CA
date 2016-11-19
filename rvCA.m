
%%

global cells cells_old n x y rulecurr async
rind=1;
k=rind;
randrule=1;
randrule=0;
name='nerd';
loadrule
getrule;

n=200;
x=2:n+1;
y=2:n+1;
cells=zeros(n+2,n+2);
async=0;
siz=size(cells);
[x1,x2]=ndgrid(x,y);
xyid=sub2ind(siz,x1,x2);
p0=0.5;

%%
cells=rand(n+2,n+2)*5-2.5;
%cells=zeros(n+2,n+2);
% cells(10:15,10:15)=rand(6,6)/2;
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
upf=@(a,b) mod(a,8)+mod(b,1);
fc='Sinput(xyid)=arrayfun(upf,Sinput(xyid),cells(xyid));';
fc='Sinput(xyid)=mod(mod(Sinput(xyid),8)./px+mod(cells(xyid),1)./py,1);';
% fc='Sinput(xyid)=(mod(Sinput(xyid),8)./px+mod(cells(xyid),1)./py)./(8./px+1./py);';
fc='Sinput(xyid)=(mod(Sinput(xyid),8)./px+mod(cells(xyid),1)./py);';
updateFcn=[fc];
% f2='Sinput(xyid);'

% px=px(128,1);
% py=py(1,96);

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
set(fi.Parent,'XTickLabels',cellstr(num2str(px(ax.XTick,1),3)));
set(fi.Parent,'YTickLabels',cellstr(num2str(py(1,ax.YTick)',3)));

while stepnum<stepmax
    %%
cold=cells(xyid);
cells=floor(div*cells)/div;
Sinput=conv2(cells,[1 1 1; 1 0 1; 1 1 1],'same');
% Sinput=(floor(div*Sinput)/div);
% cells=rulecurr(Sinput);
% Sinput(xyid)=eval(f1);
% Sinput(xyid)=eval(f2);
eval(updateFcn);
% Sinput(xyid)=mod(Sinput(xyid),px)./py;
cells(xyid)=Sinput(xyid);
cells=torus(cells);
cellsT=cells(xyid);
stepnum=stepnum+1;
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
im=frame2im(getframe(fi.Parent.Parent));
[imind,cm]=rgb2ind(im,256);

if stepnum==intl;
    imwrite(imind,cm,[gifname '.gif'],'gif','Loopcount',inf,'DelayTime',0.05);
else
    imwrite(imind,cm,[gifname '.gif'],'gif','WriteMode','append','DelayTime',0.05);
end
pause(0.05)
% stepnum=stepnum-intl+1;
end

end

%%
[min(px(:)) max(px(:));
  min(py(:)) max(py(:))]'