ac=1000;
upd=@(a,mu) mu.*a.*(1-a);
a=rand(1,ac);
a=linspace(0,1,ac);
% mu=linspace(3.8281,3.829,ac);
mu=linspace(1,4,ac);
[a,mu]=ndgrid(a,mu);


intl=2;
figure(1)
fi=imagesc(a);
ax=fi.Parent;
set(ax,'XTickLabels',cellstr(num2str(mu(1,ax.XTick)')))
figure(2)

for stepnum=1:stepmax
a=upd(a,mu);
if mod(stepnum,intl)==0
set(fi,'CData',a)
ax=fi.Parent;
set(ax,'XTickLabels',cellstr(num2str(mu(1,ax.XTick)')))

% set(s2,'XData',[s2.XData mu(:)'],'YData',[s2.YData a(:)'])
% scatter(mu(:),a(:),1,'.');
hold on
drawnow
% pause
end

if mod(stepnum,200)==0;
figure(2)
hold off
s2=scatter(mu(:),a(:),1,'.');
hold on

HDa=hist(a>0.5',0:1)';
HDa=bsxfun(@rdivide,HDa,sum(HDa,2));
H=sum(-HDa.*log2(HDa),2,'omitnan');
plot(mu(1,:),H)
end
end
%%
% HDa=hist(a>0.5',0:1)';
% HDa=bsxfun(@rdivide,HDa,sum(HDa,2));
% H=sum(-HDa.*log2(HDa),2,'omitnan');
% figure(3)
% plot(mu(1,:),H)
%%
figure(3)
histogram(a(:),100)