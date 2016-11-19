stds=stdfilt(cells);
cur1=(abs((sig(py/(mod1))-0.5)*2*sum(fir(:))+px*1)<0.1*2)';
cur2=(abs(px+mod1/mod2*py)<0.1*15)';
cur=cur1|cur2;
stds=repmat(stds(xyid)',1,1,3);
overlay1 = imoverlay(stds*2.5,cur, [.3 1 .3]);
figure(4)
subplot(3,1,[1 2])
f1=imagesc([stds*100 overlay1]);
subplot(3,1,3);
figure(6)
m1=linspace(-mod1,mod1,100);
m2=linspace(-mod2,mod2,100);
[m1,m2]=ndgrid(m1,m2);
o=mod(m1,mod1)/ppx+mod(m2,mod2)/ppy;
histogram2(S_input,cold,50);
set(gca,'ZScale','log')
hold on
surface(m1,m2,10.^(o))
% mesh(m1,m2,o)
shading interp
set(gca,'ydir','reverse')
hold off

title('fitted with (sig(py/3)-0.5)*15-px*1=0')
ax=f1.Parent;
xtick=mod(ax.XTick-1,400)+1;
ytick=repmat(ax.YTick,1,2);
set(ax,'XTickLabels',cellstr(num2str(px(xtick,1),3)));
set(ax,'YTickLabels',cellstr(num2str(py(1,ytick)',3)));
saveas(gcf,sprintf('phase_ex%d_fitted_framed.jpg',ex));
% figure(5)
% imagesc(stds(xyid)',[0 3]);
% ax=gca;
% set(gca,'XTickLabels',cellstr(num2str(px(ax.XTick,1),3)));
% set(gca,'YTickLabels',cellstr(num2str(py(1,ax.YTick)',3)));
% saveas(gcf,sprintf('phase_ex%d.jpg',ex));


colormap gray

