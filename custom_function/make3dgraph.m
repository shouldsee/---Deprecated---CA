function [ ax ] = make3dgraph( adj,XYZ,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if length(varargin)>=1
    thres=varargin{1};
else
    thres=50;
end
    N=length(XYZ);
    x = XYZ(1:N,2); y = XYZ(1:N,3); z = XYZ(1:N,1);
    thres=prctile(adj(:),thres);
%     thres=0;
    [r c] = find(adj>thres);
    p = [r c]';              %'# indices
    ax=plot3(x(p), y(p), z(p), 'LineWidth',0.5, 'Color',[.4 .4 1], ...
    'Marker','o', 'MarkerSize',2, ...
    'MarkerFaceColor','g', 'MarkerEdgeColor','g');
wts=adj(adj>thres);
wts=floor(wts/max(wts)*20)+1;
cmap=colormap(gca);
for i=1:length(ax)
%     LWd=LWds(i)+20;
    set(ax(i),'color',cmap(wts(i),:));
end

end

