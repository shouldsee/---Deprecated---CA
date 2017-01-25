lag=1;
n=200;
% name='nerd';
% loadrule
rind=1;
rulecurr=ones(1,18);
brun=0;
% getrule
init_CA
figure(3)
beta0=linspace(-5,5,n+2)/1;
bmat=repmat(1./beta0',1,n+2);
b = uicontrol('Style','slider','Min',min(beta0),'Max',max(beta0),...
                'SliderStep',diff(beta0(1:3)),'Value',beta0(1),...
                'Position',[20 20 200 20],...
                'CallBack',@(hObj,eventdata) get(hObj,'Value'));
intel = uicontrol('Style','slider','Min',1,'Max',8000,...
    'SliderStep',[0.5 0.5],'Value',4,...
    'Position',[200 20 200 20],...
    'CallBack',@(hObj,eventdata) get(hObj,'Value'));
            [outstring,newpos] = textwrap(b,cellstr(sprintf('%.2f',get(b,'Value'))));

            %     set(b,'Callback',); 
cells=(single(rand(n+2,n+2)<0.8)-0.5)*2;
fi=imagesc(cells);
fir=[0 1 0;
    1 0 1;
    0 1 0];
% fir=ones(3,3);
fir(2,2)=0;
nfir=fir/sum(fir(:));
avg=conv2(cells,nfir,'same');
siz=size(cells)-[4 4];
stepnum=1;
intl=4E1;
px=repmat(linspace(0,1,n+2),n+2,1);

ax=gca;
ytick=floor(linspace(1,n+2,20));
set(ax,'YTick',ytick,'YTickLabels',(1./beta0(ytick')));
[x1,y1]=ndgrid(1:n+2,1:n+2);
chk=xor(mod(x1,2),mod(y1,2));
% cold=rand(n+2,n+2)<0.5;
% coold=rand(n+2,n+2)<0.5;
while true;
    if mod(stepnum,intl)==0;
%         set(fi,'CData',(stdfilt(cells)));
        set(fi,'CData',(cells));
        drawnow
        brun=1/get(b,'Value');
        intl=ceil(get(intel,'Value'));
        fprintf('beta=%.2f, stepnum=%d, energy=%.3d, mag=%.3d \n',brun,stepnum,mean(avg(:).*cells(:)),mean(cells(:)));
    end
    id=randi([1 prod(siz)],[1,2]);
    [xi,yi]=ind2sub(siz,id);
    xi=xi+2;
    yi=yi+2;
    id=sub2ind([n+2,n+2],xi,yi);
    %% sequential udt
    if sum(cells(id))==0;
%    p=exp(mean(-brun.*2.*cells(id).*avg(id)));
% %    p=exp(-bmat(xi,yi)*2*cells(xi,yi)*avg(xi,yi));
% %    udt=p>1;
%    p=min(1,p);
%    udt=rand(1,1)<p;
    udt=mean(-brun.*2.*cells(id).*avg(id))>0;
%    udt=p==1;
%     udt=avg(xi,yi)==0;
   if udt;
       cells(id)=-cells(id);
       avg(xi-1:xi+1,yi-1:yi+1)=conv2(cells(xi-2:xi+2,yi-2:yi+2),nfir,'valid');
   end
   
    end
    %% sync udt
% %     pm=exp(-brun*2*cells*avg);
%     pm=exp(-bmat.*2.*cells.*avg);
%     udt=(rand(n+2,n+2)<pm);
% %     udt=pm>1;
%     cells=(udt.*(rand(n+2,n+2)<px)-0.5)*-2.*cells;
% % cells=(udt-0.5)*-2.*cells;
% %    
    %% Q2R CA
%     avg=conv2(cells,nfir,'same');
%     rev=((avg==0 & chk)-0.5)*2;
%     cells=cells.*-1.*rev;
%     chk=1-chk;
%     cells=(((avg==0)-0.5).*2.*coold);
    
    %%
    stepnum=stepnum+1;
end