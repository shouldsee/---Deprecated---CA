stds=stdfilt(cells);
cur1=(abs((sig(py/(mod1))-0.5)*2*sum(fir(:))+px*1)<0.1*2)';
cur2=(abs(px/(mod1-mod2/2)+py/mod2)<0.1*1)';
cur=cur1|cur2;
stds=repmat(stds(xyid)',1,1,3);
overlay1 = imoverlay(stds*2.5,cur, [.3 1 .3]);

figure(4)
subplot(3,1,[1 2])
f1=imagesc([stds*100 overlay1]);
subplot(3,1,3);
%%
figure(6)
clf
hh=histogram2(cold(:),cold(:),50);
hold on
mh=mesh([],[],[],[]);
hold on
idm=mesh([],[],[],[]);
shading interp
view(-40,40)
xlabel('S input')
ylabel('cells')
colormap default
xlim([-1 1])
ylim([-1 1])
title('the map function, the identity plane, and the distribution')
b = uicontrol('Style','slider','Min',-1,'Max',1,...
                'SliderStep',[0.01 0.01],'Value',0.47,...
                'Position',[20 20 200 20],...
                'CallBack',@(hObj,eventdata) get(hObj,'Value'));
number = uicontrol('style','text', ...
    'string','1', ...
   'fontsize',12, ...
   'position',[20,40,50,20],...
   'string','0.5');
            
%%
while true
   wrap1=@(cells,Sinput) mod(Sinput+bnow,1)./ppx+mod(cells,1)/ppy;

braid
% set(hh,'Data',[S_input cold(:)]);
drawnow
pause(0.2)

end
%%

% title('fitted with (sig(py/3)-0.5)*15-px*1=0')
ax=f1.Parent;
xtick=mod(ax.XTick-1,n)+1;
ytick=repmat(ax.YTick,1,2);
set(ax,'XTickLabels',cellstr(num2str(px(xtick,1),3)));
set(ax,'YTickLabels',cellstr(num2str(py(1,ytick)',3)));
saveas(gcf,sprintf('../data/phase_ex%d_fitted_framed.jpg',ex));
% figure(5)
% imagesc(stds(xyid)',[0 3]);
% ax=gca;
% set(gca,'XTickLabels',cellstr(num2str(px(ax.XTick,1),3)));
% set(gca,'YTickLabels',cellstr(num2str(py(1,ax.YTick)',3)));
% saveas(gcf,sprintf('phase_ex%d.jpg',ex));

% function updtplot(hObject,event,x,hplot)
% n = get(hObject,'Value');
% set(hplot,'ydata',x.^n);
% drawnow;
% end
% fprintf('%d\t%d\t%d\t%.4f\n',mod1,mod2,sum(fir(:)),ppx/ppy);
fprintf('%.f\t\t\t%.f\t\t\t%d\t\t\t%.4f\n',[mod1 mod2 sum(fir(:)) ppx/ppy]);
% [mod1 mod2 sum(fir(:)) ppx/ppy]
%%


% colormap(c.map);
