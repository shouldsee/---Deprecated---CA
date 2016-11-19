% load('visual_complex','S_poplist','rulename')
% S_poplist_co=S_poplist;
% rulename_co=rulename;
% load('visual_chaos','S_poplist','rulename')
% S_poplist=[S_poplist_co;S_poplist];
% rulename=[rulename_co;rulename];
% ruletemp=importrule(rulename);
% save visual_both
%%
update1=0;
update2=1;
pause1=0;
pause2=0;
store=1;
% kstep=6;
% csvread('nerd.csv')
% T=table2cell(readtable('nerd.csv','Delimiter',' '));
% rulename=T(:,1);
% id=cellfun(@(t)~strcmp(t,''),rulename);
% rulename=rulename(id);
% ruletemp=importrule(rulename);
% load(name,'rulename');

if~exist('galleryname')
galleryname='gallery/'
end
if ~exist('batchname')
batchname='temp/';
end
if ~exist('vmethod')
    vmethod='rp_visual6'
end

galleryname=['gallery/' vmethod '/'];
dirname=[galleryname batchname];
if ~exist(dirname,'dir')
    mkdir(dirname);
end
% k0=1;
eval(vmethod);
%%
% figure(10)
% v1s_bak=v1s;
% v2s_bak=v2s;
% % subplot(1,2,1)
% data=v1s_bak;
% xs=data(:,1);
% ys=data(:,2);
% c=[ones(1,40) 10*ones(1,40)];
% dict_legend={'complex','chaotic'};
% leg=dict_legend{c};
% gscatter(xs,ys,c)
% 
% text(xs,ys,num2str((1:length(xs))'))
% title('soups scattered on space of complexity')
% legend('complex','chaotic','Location','Best');
% xlabel('mean(stat entropy)')
% ylabel('variance(stat entropy) expained by trend')



% 
% 
% subplot(1,2,2)
% data=v2s_bak;
% xs=data(:,1);
% ys=data(:,2);
% scatter(xs,ys)
% text(xs,ys,num2str((1:length(xs))'))

