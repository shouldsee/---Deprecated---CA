
global cells cells_old n x y
n=30;
x=2:n+1;
y=2:n+1;
cells=zeros(n+2,n+2);
% cells(2:pop+1,2:pop+1)=single(rand(pop,pop)<p0);
soup=single(rand(n,n)<p0);
% soup=glider;
cells(x,y)=soup;


step=1;
dfx=randi([2 n+1]);
dfy=randi([2 n+1]);
dfx0=dfx;
dfy0=dfy;
df_old=zeros(lag+1,2);
time=0;
hop=0;
dso=[];
movx=zeros(pop,pop);
movy=zeros(pop,pop);
ds=zeros(pop,pop);
%%
k=k0;
rind=floor((k-1)/soupdensi)+1;    
getrule
% soup=1;
%%
% cells=update(cells,rulecurr);
%%
% subplot(2,2,3)
% f3s=scatter([],[],3 ,'x');
clear s
s.dfx=1;
s.dfy=1;
s.dfx0=1;
s.dfy0=1;
s.ddfx=0;
s.ddfy=0;
s.ds=0;
s.pt=0;
% s.time=[];
s.ind=1;

particle=repmat(s,n,n,lag+1,2);
siz=size(particle);
xs_old=zeros(n,n,lag+1,lag+1);
ys_old=zeros(n,n,lag+1,lag+1);
%%
subplot(2,2,1)
f1s=scatter([],[],3,'x');
xlim([0 20])
ylim([0 20])
subplot(2,2,3)
f3h=histogram2([],[],(-lag-1:lag+2)-0.5,(-lag-1:2*lag+2)-0.5);
set(gca,'ZScale','log')
line(1:2*lag+2,1:2*lag+2,ones(2*lag+2,1));
line(1:2*lag+2,(1:2*lag+2).^0.5,ones(2*lag+2,1));

middle=x(ceil(numel(x)/2));
view(60,45)
xlabel('xaxis')
ylabel('yaxis')


subplot(2,2,3)
f3i=imagesc([]);
caxis([0 10])
subplot(2,2,4)
f4i=imagesc([]);
caxis([0 10])

% cells=back;
% f4s=scatter([],[],3,'x');
while step<1000
% cellsdftd=cells;
cells_old=cells;
cellsdf=cells_old;
cells=update(cells);
step=step+1;    

for p1=middle:middle+pop-1
    for p2=middle:middle+pop-1
        for pt=size(particle,3):-1:1
            for ps=1:size(particle,4)
                p=particle(p1,p2,pt,ps);
                if pt==lag+1
                    movx(p1,p2)=(p.dfx-p.dfx0);
                    movy(p1,p2)=(p.dfy-p.dfy0);
                    ds(p1,p2)=p.ds;
                    p=particle(p1,p2,pt-1,ps);
                elseif pt==1
                    p.dfx0=p1+1;
                    p.dfy0=p2+1;
                    p.dfx=p.dfx0;
                    p.dfy=p.dfy0;
                    p.ds=0;
                else
                    p=particle(p1,p2,pt-1,ps);
                end
                cellsdf(p.dfx,p.dfy)=cellsdf(p.dfx,p.dfy)+5;

                p=update_p(p);
                
                    particle(p1,p2,pt,ps)=p;
            end
        end
    end
end

subplot(2,2,2)
imagesc(cellsdf(x,y))
% xs=sum(abs([movx(:),movy(:)]),2);
% ys=ds(:);xs=reshape([particle.pt],n,n,lag+1);


xs0=cat(5,reshape([particle.ddfx],siz),reshape([particle.ddfy],siz));
xs=squeeze(sum(abs(diff(xs0,[],4)),5));
xs_old=cat(4,xs,xs_old(:,:,:,1:lag));

ys0=cat(5,reshape([particle.ddfx],siz),reshape([particle.ddfy],siz));
ys=squeeze(sum(abs(mean(ys0,4)),5));
ys_old=cat(4,ys,ys_old(:,:,:,1:lag));


% xs=convn(xs,[1 1 1;1 -9 1;1 1 1]/9,'same');
% xs=xs(:,:,lag+1,:);
% xs=sum(abs(xs),4);
% ys=convn(ys,[1 1 1;1 1 1;1 1 1]/9,'same');
% ys=ys(:,:,lag+1,:);
% ys=sum(abs(ys),4);


% % set(f1s)


f3im=squeeze(xs_old(middle,middle,:,:));
f4im=squeeze(ys_old(middle,middle,:,:));

set(f1s,'XData',f3im(:),'YData',f4im(:));
set(f3i,'CData',f3im);

set(f4i,'CData',cells(x,y));

% imagesc(ys)
% caxis([0 10])

% xs=reshape([particle.pt],n,n,lag+1);
% xs=sum(sum(xs,1),2);
% ys=reshape(abs([particle.dfx]-[particle.dfx0])+abs([particle.dfy]-[particle.dfy0]),n,n,lag+1);
% ys=sum(sum(ys,1),2);
% set(f4s,'XData',xs(:),'YData',ys(:));
% set(f3h,'Data',[xs(:),ys(:)]);
if mod(step,2)==0
drawnow
end
%             if step==1||time ==10
%             end
    
    % dfds=xor(update(cells,rulecurr),update(cellsdftd,rulecurr));

%     time=time+1;    
end
% subplot(2,2,1)
% imagesc(dfds(x,y))

% subplot(2,2,2)
% imagesc(cellsdf(x,y))
% % xs=[f3s.XData df_old(2,1)];
% % ys=
% set(f4s,'XData',xs,'YData',[f3s.YData df_old(2,2)]);
% 
% drawnow
% if isempty(ind_dfds)||time>5
% 
%   if time>5
%       time=0;
%       ds=[dfx-dfx0,dfy-dfy0];
%         dso=[dso;ds];
%         dfx0=dfx;
%         dfy0=dfy;
%   end
% dfx=randi([2 n+1]);
% dfy=randi([2 n+1]);
% hop=1;
% else
% ind=datasample(ind_dfds,1);
% [dfx,dfy]=ind2sub([n n],ind);
% dfx=dfx+1;
% dfy=dfy+1;
% end
% 
% df=[dfx,dfy];
% df_old=[df;df_old(1:lag,:)];
% ddf_old=diff(df_old,[],1);

% if hop==1;
% ddf_old(:,1)=0;    
% hop=0;
% end




%%

