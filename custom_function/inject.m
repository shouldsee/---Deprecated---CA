function[gx,gy,rg]=inject(wrap1)
inj=linspace(-5,5,1000);
[gx gy]=ndgrid(inj,inj);
gz=wrap1(gx(:),gy(:));
MIN=min(gz);
MAX=max(gz);
rg=[MIN,MAX];
rs=linspace(MIN,MAX,50);
[gx,gy]=ndgrid(rs,rs);
end