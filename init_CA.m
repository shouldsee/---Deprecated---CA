randrule=0;
rind=1;
% n=200;
x=2:n+1;
y=2:n+1;
[x1,y1]=ndgrid(x,y);
xyid=sub2ind([n+2,n+2],x1,y1);
rulenum0=1;
intl=2;
dm=0;
stepmax=500;

preallocate;