global cells cells_old n x y rulecurr async
rind=1;
k=rind;
randrule=1;
randrule=0;
name='nerd';
loadrule
getrule;

n=16;
x=2:n+1;
y=2:n+1;
cells=zeros(n+2,n+2);
async=0;
siz=size(cells);
[x1,x2]=ndgrid(x,y);
xyid=sub2ind(siz,x1,x2);
p0=0.5;

cooldown=10;
soupmax=300;
rulemax=1E100;


% frames=zeros([numel(xyid),soupmax]);
frame1=cell(1,soupmax);
frame2=frame1;
%%
soupnum=1;
while soupnum<=soupmax
%     soupind=soup;
    cells(xyid)=rand(size(xyid))<p0;
    for i=1:cooldown
        cells=update(cells);
    end
    cells=torus(cells);
    S_inputT=conv2(cells,[1 1 1;1 9 1;1 1 1],'same');
    frame1{soupnum}=reshape(S_inputT(xyid),1,[]);
    cells=update(cells);
    frame2{soupnum}=reshape(cells(xyid),1,[]);

%     frames(:,soupnum)=reshape(cells(xyid),1,[]);
    soupnum=soupnum+1;
end