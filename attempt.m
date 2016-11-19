% scatter3(so(:,1),so(:,2),so(:,3),3,'x');
work=1;
%%
so=[];
so.stats=[];
so.rulelist=[];
so.label=[];
%%
% so=[];
load visual_test;
lag=2;
label=1;
scanthrough
%%
load ambig_und;
label=2;
lag=2;
scanthrough
%%
load visual_chaos;
label=3;
lag=2;
scanthrough
%%
load life
lag=2;
label=4;
scanthrough
%%
load weirdo
lag=2;
label=5;
scanthrough
%%
load ambig2
lag=2;
label=6;
scanthrough
%%

load visual_complex
lag=2;
label=7;
scanthrough

%%
%%

urule=unique(so.rulelist);
so.MapObj=containers.Map(urule,1:length(urule));
so.group=zeros(1,length(so.rulelist));
for i=1:length(so.rulelist)
    so.group(i)=so.MapObj(so.rulelist{i});
end
save('so_SI','so')


work=0;

%%
% so_ambig=[so,repmat(1,size(so,1),1)];
% so_ambig=so;
% %%
% list={'so_chaos','so_test','so_complex','so_ambig'};
% list={'so_test','so_ambig'};
% %%
% clf
% for i=1:length(list)
%     eval(sprintf('so=%s;',list{i}));
%     scatter3(so(:,1),so(:,2),so(:,3),3,'x');
%     hold on 
% end
% hold off
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,rulelist})
