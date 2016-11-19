
uni=100;
%%
S_pop=rand(n,n,uni)<p0;
h=histogram(pdist(reshape(S_pop,n^2,[]),'hamm'),0:0.005:1);
ylim([0 4E5])
xlim([0 1])
cellsM=zeros(n+2,n+2,uni);
i=1;

%%
while true
cellsM(x,y,:)=S_pop;
cellsM(1,:,:)=cellsM(x(end),:,:);
cellsM(end,:,:)=cellsM(x(1),:,:);
cellsM(:,1,:)=cellsM(:,y(end),:);
cellsM(:,end,:)=cellsM(:,y(1),:);
S_inputT=convn(single(cellsM),FIR.S_input,'same');
S_input=S_inputT(x,y,:);
S_pop=rulecurr(S_input+1);
d=pdist(reshape(S_pop,n^2,[]),'hamm');
set(h,'Data',d);
phase(i)=mean(d);
i=i+1;
drawnow
% pause(0.2)
end