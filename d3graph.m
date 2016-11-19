%# sample adjacency matrix and 3D coordinates of points
% N = 30;                                      %# number of vertices
% [adj,XYZ] = bucky;
% adj = full(adj); adj = adj(1:N,1:N);

N=length(V_input);
% adj=dmat_V_input;
adj=(MD_Tinputo+MD_Tinputo')/2;
XYZ=V_input(:,1:3);
labels = cellstr( num2str((1:N)','%02d') );  %'# nodes labels
x = XYZ(1:N,2); y = XYZ(1:N,3); z = XYZ(1:N,1);
%# another sample data
%#x = rand(N,1);          %# x-coords of vertices
%#y = rand(N,1);          %# y-coords
%#z = rand(N,1);          %# z-coords
%#adj = rand(N,N)>0.7;    %# adjacency matrix

%# plot graph in 3D
[r c] = find(adj>prctile(adj(:),50));

%%
data=reshape(V_inputo,[],18,4);


%%
% subplot(3,2,1)
clf
axis([0 4 0 4 0 0.5 ]),box on
view(50,37.5);
colorbar
axs=[];
inds=[];
lag=4+1;
for j=1:length(data);
    delete(axs(inds==lag))
    inds(inds==lag)=[];
    XYZ=squeeze(data(j,:,1:3));
    adj=reshape(D_Tinputo(j,:),18,[])+20;
    ax=make3dgraph(adj,XYZ,90);
    axs=[axs;ax];
    inds=[inds+1;ones(length(axs),1)];
    axis([0 4 0 4 0 0.5 ])
%     view(50,37.5);
%     thres=prctile(adj(:),0);
hold on
% 
drawnow
% pause
end
% set(ax)
% text(x, y, z, labels, ...
%     'EdgeColor','g', 'BackgroundColor',[.7 1 .7], ...
%     'VerticalAlignment','bottom', 'HorizontalAlignment','left')
%  axis vis3d, box on, view(3)

%%
adj=(MD_Tinputo+MD_Tinputo')/2;
G=graph(adj);
LWidth=(G.Edges.Weight+20);
LWidth=LWidth/max(LWidth)*5;
p=plot(G,'LineWidth',LWidth);
xlabel('x'), ylabel('y'), zlabel('z')
p.XData=V_input(:,2);
p.YData=V_input(:,3);