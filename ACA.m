global n x y async cells rulecurr a b




%%
%% implement rules along a gardient of parameter.
%% input:
% iter:rule that map cells to next cells
% parameter to be varied in the rule
% the gradient of the given parameter
% 
% iter=@(cells,gd); 
% 
% iter=iterbase(cells,async); 
% cells=iter(cells)


figure(2)
n=400;
intl=40;
stepmax=2^18;

% async=0;
async=0.35;
a=n+2;
b=n+2;
a=1600;
b=40;
% b=40+2;
% % % % % % 
% stepmax=200;
async=1;
cells=single(rand(a,b)<p0);
cells=gpuArray(single(rand(a,b)<p0));

gda1=0.98;
gda2=0.99;
gdb1=0.98;
gdb2=0.99;
gda1=0.01;
gda2=0.99;

[gda,gdb]=ndgrid(linspace(gda1,gda2,a),linspace(gdb1,gdb2,b));
J=gda;
udrate=gdb;
% subplot(2,1,1)
% fi=imagesc(cells(xyid)');
% subplot(2,1,2)
% fi2=imagesc(gdb');
% % axis auto
fi=imagesc(cells(xyid)');
ylabel('update rate, gda')
xlabel('coupling energy, gdb')


x=2:a-1;
y=2:b-1;
[x1,x2]=ndgrid(x,y);
% p0=0.991;
p0=0.5;
m=5;
xyid=sub2ind(size(cells),x1,x2);




rtp=str2rule('B1678/S236');
rulecurr=rtp;
% rulecurr=zeros(1,18);
% rulecurr([4,12,13])=1;
% rindmax=max(size(rulename));
rindmax=1000;
fig=gcf;
% fig.Position=[-500 -500 800 200];
% figure(3)
% l=plot([]);
% rind0=1;
rindmax=2^17;
for rind=rind0:rindmax
rindbk=rind;
%     rind=rss(rind);
    rind=id(rind);
    %     for rind=1:2^17;
%%
tic
    % getrule;
    rulecurr=str2num(dec2base(rind-1,2,18)')';
-0.5,1.5
% rind;
% fig=gcf;
if fig.Number~=2;
    figure(2)
end
if ~ishandle(fi)
    fi=imagesc(cells(xyid)');
end
% cells=single(rand(a,b)<p0);

% rname=rulename{rind};
rb=find(rulecurr(1:9))-1;
rs=find(rulecurr(10:18))-1;
B=num2str(rb')';
S=num2str(rs')';
rname=['B' B '/S' S];
suptitle(sprintf('%s',rname));

for step=1:stepmax
%%
    %     cells_old=cells;
        cells=iterbasebase(1,cells,@(cells) cell2transi_ising(cells,J),udrate); 
%         cells=iter(cells);
        S_input=conv2(cells,FIR.S_input,'same');
        set(fi,'CData',gather(S_input(xyid)'));
        caxis(fi.Parent,[0 17])
        
%         set(fi,'CData',cells(xyid)')
%         caxis(fi.Parent,[0 1])
        if mod(step,intl)==0
          drawnow
        end
%         pause(0.2)
%         step
end
k=rind;
% figure(3)
h=hist(cells(xyid)',0:1);
% h=hist(S_input(xyid)',0:17);
hm=movmean(h',m,'Endpoints','discard');
HD=bsxfun(@rdivide,hm,b*m);
H=sum(-log2(HD).*HD,2,'omitnan');
rec(rind,1:length(H))=H;

% plot(H)
% ylim([0 1])

% drawnow
% rind
% makeprint
% toc
% t(rind)=toc;
%%
end

%%
rec=zeros(rindmax,398);

%%
figure(4)
imagesc(rec(1:rind,:)');

