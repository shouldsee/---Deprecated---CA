global n x y async cells rulecurr a b

figure(2)
% n=400;
n=600;

stepmax=2^18;

% async=0;
async=0.35;

a=n+2;
b=a;
% % % % % % 
n=400;
b=40+2;
stepmax=200;
async=1;
cells=single(rand(a,b)<p0);


x=2:a-1;
y=2:b-1;
[x1,x2]=ndgrid(x,y);
% p0=0.991;
p0=0.5;
m=5;
xyid=sub2ind(size(cells),x1,x2);
% cells=zeros(n+2,n+2);
sync=linspace(-0.5,1.5,n+2);
% sync(250)
fi=imagesc(cells(xyid)');
% axis auto

rtp=importrule({'B1678/S236'});
rulecurr=[rtp{1},rtp{2}];
% rulecurr=zeros(1,18);
% rulecurr([4,12,13])=1;
rindmax=max(size(rulename));
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
    %     cells_old=cells;
        cells=update(cells);
% 
        S_input=conv2(cells,FIR.S_input,'same');
        set(fi,'CData',S_input(xyid)');
        caxis(fi.Parent,[0 17])
        
%         set(fi,'CData',cells(xyid)')
%         caxis(fi.Parent,[0 1])
        if mod(step,2)==0
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

