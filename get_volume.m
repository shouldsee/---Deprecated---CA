function[V]=get_volume(cur_path);
[p,xi]=ksdensity(cur_path','Kernel','box','bandwidth',0.125 ...
    ...,'Support',[0 -inf;2*pi inf]...
    ..., 'PlotFcn','contour'    ...
    );
siz=size(xi);
msh=reshape(xi,sqrt(siz(1)),sqrt(siz(1)),siz(2));
dV=prod(squeeze(msh(2,2,:)-msh(1,1,:)));
V=sum(p(:)>0)*dV;
end