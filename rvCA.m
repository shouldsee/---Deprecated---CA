
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
if ~exist ('crg','var')
    crg=[-1 1]*1;
end
if ~exist('ex','var')
ex=15;
end

global cells cells_old n x y rulecurr async px py
rind=1;
k=rind;
% randrule=1;
% randrule=0;
name='nerd';
% loadrule
if randrule
getrule;
end

n=500;
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
cells=rand(n+2,n+2)*ex;
% cells=zeros(n+2,n+2);
cells(10:15,10:15)=rand(6,6)/2;
end
rsoup=rand(n+2,n+2)*ex;
cells=rsoup;
speed=rand(n+2,n+2)*ex;
%
cells=torus(cells);
rulecurr=@(a) 1-(mod(100*a,100)/100);
% rulecurr=@(a) (mod(a,4))/4;
rulecurr=@(a) tanh(a/10).*sin(a);

figure(1)
fi=imagesc(cells(xyid)',crg);
figure(2)
h=histogram(cells(xyid));
figure(3)
h2=histogram2(cells(xyid),cells(xyid),100);
% edge=-1:0.025:1;
edge=linspace(crg(1),crg(2),100);
h2=histogram2(cells(xyid),cells(xyid),edge,edge);
view(30,30);
zlim([10,2000]);
set(gca,'ZScale','log')

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



% f2='Sinput(xyid);'
upf=@(a,b) mod(a,8)+mod(b,1);
qmap=@(a,b) b.*a.*(1-a);
fc='Sinput(xyid)=arrayfun(upf,Sinput(xyid),cells(xyid));';
fc='Sinput(xyid)=mod(mod(Sinput(xyid),8)./px+mod(cells(xyid),1)./py,1);';
% fc='Sinput(xyid)=(mod(Sinput(xyid),8)./px+mod(cells(xyid),1)./py)./(8./px+1./py);';
fc='Sinput(xyid)=(mod(Sinput(xyid),mod1)./px+mod(cells(xyid),mod2)./py);';
% syms a b c d
% wrap_syms=mod(a,mod1)/c+mod(b,mod2)/d;
wrap=@(cells,Sinput) mod(Sinput,mod1)./px+mod(cells,mod2)./py;
% wrap1=@(cells,Sinput) mod(Sinput,px)+mod(cells,px);
% wrap1=@(cells,Sinput) mod(Sinput,1)./px+mod(cells,1)./py;
ppx=-1.2;
% ppy=3.3;
ppy=5.86
wrap1=@(cells,Sinput) mod(Sinput,1)./ppx+mod(cells,1)/ppy;

% T = convmtx2(H,m,n) 
% reshape(T*cells(:),size(H)+[m n]-1); 
fc=['Sinput=conv2(cells,nfir,''same'');\n',...
    'Sinput(xyid)=wrap1(cells(xyid),Sinput(xyid));\n',...
    'cells(xyid)=Sinput(xyid);',...
    'cells=torus(cells);'],...
%     'lyap='];
% fc=['cellsT=cells(xyid);',...
% 'Sinput=reshape(TT*cellsT(:),n,n);',...
% 'cells(xyid)=wrap(Sinput,cells(xyid));',...
% 'jcb=(ddpx+ddpy)*TT*jcb;'];

% fc=['cells(xyid)=qmap(cells(xyid),1.*px);\n',...
%     'cells=torus(cells);\n',...
%     'Sinput=conv2(cells,nfir,''same'');\n',...
%     'S_input=Sinput(xyid);\n',...
%     'Sinput(xyid)=(Sinput(xyid)-cells(xyid)).*py+(1-1.*py).*cells(xyid);\n',...
%     'Sinput=max(min(Sinput,1),0);\n',...
%     ];
% fc=['cells=torus(cells);\n',...
%     'speed=torus(speed);\n',...
%     'speed(xyid)=speed(xyid)-(conv2(cells(xyid),nfir,''same'')-cells(xyid))./px;\n',...
%     'Sinput(xyid)=cells(xyid)+speed(xyid);\n',...
% %     'Sinput(xyid)=(Sinput(xyid)-cells(xyid)).*py+(1-1.*py).*cells(xyid);\n',...
% %     'Sinput=max(min(Sinput,1),0);\n',...
%     ];
fc=sprintf(fc);
%     'Sinput(xyid)'
updateFcn=[fc];
if ~exist('mod1','var')
mod1=4;
end
if ~exist('mod2','var')
    mod2=1;
end
% px=3.6*ones(size(xyid));
% py=0.7*ones(size(xyid));

% py=py/mod1*mod2;
% py=-px./((mod1-mod2/2)).*(mod2)+py./ex.*px./ex/0.2;
% py=-px./(mod1-1).*(mod2);
% rect=[30,224;
%       30,227]; %bb-mosa interface
% zoomin([30,224;30,227])
% zoomin([15,207;32,221])
% zoomin([90,283;96,290])
% zoomin([110,310;126,326])
% zoomin([86,326;92,332])
% zoomin([81,382;92,392])
% zoomin([201,201;400,400])
% rect=[100,180;
%       103,215]; %bb-mosa interface
% rect=[100,210;
%       103,215]; %bb

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
% % rect=[300,140;329,155];
% rect=[163,60;164,60];


% rect=[80,208;
%       128,240];
% rect=[90,180;
%       128,380];
% rect=[80,240;
%       80,240];

% px=repmat(linspace(px(rect(1,1),1),px(rect(2,1),1),n)',1,n);
% py=repmat(linspace(py(1,rect(1,2)),py(1,rect(2,2)),n),n,1);

rect=[300,140;329,155];
rect=[163,60;164,60];

% rect=[193,209;194,210];
% rect=[200,96;200,96];
% rect=[51,64;126,121];

% rect=[159,283;
%       159,283];

% % rect=[100,50;
% %       103,180];
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

ax=fi.Parent;
fir=[1 1 1; 1 0 1; 1 1 1];
% fir=[1 1 1; 1 1 1; 1 1 1];
fir=[0 1 0; 1 1 1; 0 1 0];
fir(2,2)=0;
% fir=[0 1 0; 1 1 1; 1 1 0];
% fir=ones(3,3)*3;
% fir=rand(5,5);

nfir=fir/sum(fir(:));
% TT=linconv(nfir,cells(xyid));
% ddpx=sparse(diag(1./px(:)));
% ddpy=sparse(diag(1./py(:)));
% TT=sparse(TT);
% jcb=((ddpx+ddpy)*TT);

% frame=speye([n,n,n,n]);
% % ddpx=sqrt(bsxfun(@times,ddpx,shiftdim(ddpx,-2)));
% ddpy=sqrt(bsxfun(@times,ddpy,shiftdim(ddpy,-2)));
% ddpx=sparse(reshape(ddpx,n^2,n^2));
% ddpy=sparse(reshape(ddpy,n^2,n^2));

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
figure(1)
% cells=rand(n+2,n+2)+0.1;
for stepnum=1:stepmax
    %%
    pr=rand(n,n);
% cells=floor(div*cells)/div;
cold=gather(cells(xyid));

% S_inputold=S_input;

% Sinput=(floor(div*Sinput)/div);
eval(updateFcn);
% cells=rulecurr(Sinput);
% Sinput(xyid)=eval(f1);
% Sinput(xyid)=eval(f2);

% % Sinput(xyid)=mod(Sinput(xyid),px)./py;
% cells(xyid)=Sinput(xyid);
% cells=torus(cells);
% cellsT=gather(cells(xyid));
% %stepnum=stepnum+1;
if mod(stepnum,intl)==0
set(fi,'CData',gather(cells(xyid)'))
% set(fi,'CData',stdfilt(cells(xyid))');
set(h,'Data',gather(cells(xyid)))
% set(h2,'Data',[cold(:),cellsT(:)])
% set(h2,'Data',gather([S_input(:),cold(:)]));
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
pause(0.01)
% stepnum=stepnum-intl+1;
end
mvs(stepnum)=gather(mv);
end

%%
fitted
% figure(6)
% histogram(cells(xyid),linspace(0,0.4,100));
% figure(7);
% plot(linspace(0,1,length(mvs)),abs(fft(mvs)))
%%

(px(:)); min(px(:))
[min(px(:)) max(px(:));
  min(py(:)) max(py(:))]'