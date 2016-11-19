cells=single(rand(n+2,n+2)<p0);

step=1;
dfx=randi([2 n+1]);
dfy=randi([2 n+1]);
dfx0=dfx;
dfy0=dfy;
df_old=zeros(lag+1,2);
time=0;
hop=0;
dso=[];
%%

%%
cells=update(cells,rulecurr);
%%
subplot(2,2,3)
f3s=scatter([],[],3 ,'x');
clear s
s.dfx=[];
s.dfy=[];
s.dfx0=[];
s.dfy0=[];
% s.time=[];
s.ind=[];
particle=repmat(s,n,n,lag+1);
%%
while step<10000
% cellsdftd=cells;
cells_old=cells;
cellsdf=cells_old;

cells=update(cells,rulecurr);

step=step+1;
    
for p1=1:size(particle,1)
    for p2=1:size(particle,2)
        for pt=size(particle,3):-1:1
            p=particle(p1,p2,pt);
            if pt==lag+1
                movx(p1,p2)=p.dfx-p.dfx0;
                movy(p1,p2)=p.dfy-p.dfy0;
            elseif pt==1
                p.dfx0=p1+1;
                p.dfy0=p2+1;
                p.dfx=p.dfx0;
                p.dfy=p.dfy0;
            else
                p=particle(p1,p2,pt-1);
            end
            
            
            particle(p1,p2,pt)=p;
        end
    end
end
%             if step==1||time ==10
%             end
    
    cellsdf(dfx,dfy)=cells_old(dfx,dfy)+10;
    cellsdftd=cells_old;
    cellsdftd(dfx,dfy)=1-cells_old(dfx,dfy);
    cellsdftd=update(cellsdftd,rulecurr);
    dfds=xor(cells,cellsdftd);
    % dfds=xor(update(cells,rulecurr),update(cellsdftd,rulecurr));
    
    ind_dfds=find(dfds(x,y));
    if isempty(ind_dfds)
        dfx=dfx;
        dfy=dfy;
    else
        ind=datasample(ind_dfds,1);
        [dfx,dfy]=ind2sub([n n],ind);
        dfx=dfx+1;
        dfy=dfy+1;
    end
%     time=time+((norm(df-df_old(end,:),2))<1);
    cellsdf(dfx,dfy)=cells(dfx,dfy)+5;
    time=time+1;

    
    
end

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

subplot(2,2,1)
imagesc(dfds(x,y))
subplot(2,2,2)
imagesc(cellsdf(x,y)+2*dfds(x,y))
set(f3s,'XData',[f3s.XData df_old(2,1)],'YData',[f3s.YData df_old(2,2)]);
drawnow
end


%%

