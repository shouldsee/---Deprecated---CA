
depth=4;
cells=single(rand(n+2,n+2)<0.5);
cellsM=repmat(cells,1,1,n+2,n+2);
cellsM_old=repmat(cellsM,1,1,1,1,depth+1);
cells_old=repmat(cells,1,1,1,1,1,depth+1);
x=2:n+1;
y=2:n+1;
%%
[s1,s2]=ndgrid(x,y);
mind=sub2ind(size(cellsM),s1,s2,s1,s2);
cmap=colormap('colorcube');
colormap default

%%
for k=k0:kstep:kmax
        rind=ceil(k/soupmax);
    suptitle(['MISC' ' rulenumber='  num2str(rind) ' ' rulename{rind}])
    r=ruletemp(rind,:);
    rulecurr=[r{1},r{2}];
% rulecurr=ruletemp
i=0;
cells=single(rand(n+2,n+2)<0.5);
phase=[];
% stepmax=depth+1
while i<20;

    S_inputT=conv2(cells,FIR.S_input,'same');
    S_input=S_inputT(x,y);
    cells(x,y)=rulecurr(S_input+1);                                    
    cells(1,:)=cells(x(end),:);
    cells(end,:)=cells(x(1),:);
    cells(:,1)=cells(:,y(end));
    cells(:,end)=cells(:,y(1));

    
    S_inputMT=convn(single(cellsM_old),FIR.S_input,'same');
    S_inputM=S_inputMT(x,y,:,:,:);
    cellsM_old=rulecurr(S_inputMT+1);
%     LS_popM=squeeze(reshape(cellsM_old(x,y,:),n^2,1,[]))';
    
cellsM=repmat(cells,1,1,n+2,n+2);
cellsM(mind)=1-cellsM(mind);
cells_old=cat(6,cells,cells_old(:,:,:,:,:,1:depth));
cellsM_old=cat(5,cellsM,cellsM_old(:,:,:,:,1:depth));


    cellsM_old(1,:,:,:,:)=cellsM_old(x(end),:,:,:,:);
    cellsM_old(end,:,:,:,:)=cellsM_old(x(1),:,:,:,:);
    cellsM_old(:,1,:,:,:)=cellsM_old(:,y(end),:,:,:);
    cellsM_old(:,end,:,:,:)=cellsM_old(:,y(1),:,:,:);

d=bsxfun(@ne,cells_old,cellsM_old);
dr=reshape(d(x,y,x,y,:,:),n^2,1,n,n,depth+1,[]);
dm=squeeze(sum(dr,1));
dv=1:depth+1;
dm=bsxfun(@rdivide,dm,shiftdim((2*(dv-1)+1).^2,-1));
vid=(dm(:,:,depth+1,1));
set(f1,'CData',cells_old(x,y,:,:,:,2));
set(f2,'CData',dm(:,:,depth+1,1));
% set(h,'Data',vid(:));
drawnow
i=i+1;
phase(i,:)=[sum(vid(:)) std(vid(:))];

end

phase(1:depth,:)=[];
phase=movmean(phase,4);
line(ax1,phase(:,1),phase(:,2),'Color',cmap(rind,:));
% line(ax1,1:size(phase,1),phase(:,1),'Color',cmap(rind,:),'DisplayName',num2str(rind));
text(ax1,phase(end,1),phase(end,2),rulename{rind},'Color',cmap(rind,:));
% line(ax2,1:size(phase,1),phase(:,2),'Color',cmap(rind,:),'DisplayName',num2str(rind));
% text(ax2,size(phase,1),phase(end,2),rulename{rind},'Color',cmap(rind,:));
% pause
end
%%