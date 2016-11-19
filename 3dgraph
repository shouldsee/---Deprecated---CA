%# sample adjacency matrix and 3D coordinates of points
% N = 30;                                      %# number of vertices
% [adj,XYZ] = bucky;
% adj = full(adj); adj = adj(1:N,1:N);
% x = XYZ(1:N,1); y = XYZ(1:N,2); z = XYZ(1:N,3);
N=length(V_inputo);
adj=dmat_V_inputo;
XYZ=V_inputo(:,1:3);
labels = cellstr( num2str((1:N)','%02d') );  %'# nodes labels

%# another sample data
%#x = rand(N,1);          %# x-coords of vertices
%#y = rand(N,1);          %# y-coords
%#z = rand(N,1);          %# z-coords
%#adj = rand(N,N)>0.7;    %# adjacency matrix

%# plot graph in 3D
[r c] = find(adj>prctile(adj(:),80));
p = [r c]';              %'# indices
plot3(x(p), y(p), z(p), 'LineWidth',2, 'Color',[.4 .4 1], ...
    'Marker','o', 'MarkerSize',10, ...
    'MarkerFaceColor','g', 'MarkerEdgeColor','g')
text(x, y, z, labels, ...
    'EdgeColor','g', 'BackgroundColor',[.7 1 .7], ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left')
axis vis3d, box on, view(3)
xlabel('x'), ylabel('y'), zlabel('z')