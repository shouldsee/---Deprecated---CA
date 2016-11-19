
% T1=readtable([fname '.xls']); %plug your stats filename into here
% T1=readtable('AICba.xls');
%%
%%
figure(10)
clf
%%
% T1=readtable('tst1.xls');

% T1=readtable('pcami_input_A.xls');
% T1=readtable('perlocate_input_C.xls');
name='flow_input_C.xls';
T1=readtable(name);
makesummary

dcm_obj = datacursormode(gcf);
r=repmat(rulelist,2,1);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,r})

hold on

%%
% T1=readtable('ambig_test3.xls');
% T1=readtable('vic_test1.xls');
% T1=readtable('cov_test1.xls');
% T1=readtable('life_test1.xls');
% T1=readtable('complex_test1.xls');
% T1=readtable('B045678_test1.xls');
% T1=readtable('vic_input_BB.xls');
% T1=readtable('percolate_visual_test.xls');
name='flow_nerd.xls';
T1=readtable(name);
makesummary

xlim([0 4])
% ylim([0 5])


% T1=readtable('pcami_weirdo.xls');
% T1=readtable('pcami_visual_test.xls');
% T1=readtable('pcami_visual_test.xls');
% makesummary

% rulelist=table2cell(T1(:,1));
%%
% T1=[T1];
% T2=[T2];
% g=[zeros(height(T1),1);1*ones(height(T2),1)];
% T1=[T1;T2];
%%
% colormap(gca,'colorcube')
% stats=table2array(T1(:,2:end-1));
% ind=find(stats(:,9)>100);
% stats=stats(ind,:);
% estep=stats(:,end);
% list1=stats(:,7);
% list2=stats(:,8);
% list2=stats(:,6)./stats(:,7);
% list3=stats(:,9);
% list3=1:length(list1);
% list4=stats(:,4);
% list5=stats(:,5);
% list6=stats(:,6);
% list7=stats(:,7);
% list8=stats(:,8);
% scatter3(list1,list2,list3,3,'x');
% plot3(list1,list2,list3);

% gscatter(list2,list3,g,colormap(gca),'x',5);
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,rulelist})
% zlim([0 0.25])

% scatter3(list1,list2,list3);
% scatter(list2,list3,3,'x');

% scatter3(log2(list1),log2(list2),log2(list3))
% hold on
% stats=table2array(T2(:,2:7));
% list1=stats(:,1);
% list2=stats(:,2);
% list3=stats(:,3);
% list4=stats(:,4);
% list5=stats(:,5);
% list6=stats(:,6);
% xlabel('mean flux(FE+BE)')
% ylabel('mean flux(FE-BE)')

% xlabel('mean[ H(input-last2) ]');
% ylabel('mean[ spatial I(input-last2) ])');
% zlabel('effective steps')
xlabel('x')
ylabel('y')
zlabel('z')
% zlabel('cov[ H(NH) , I(NH_0,NH_-_4) ]');

legend({'randomly generated-rules','curated complex rules'},'Location','best')
title('B/S rulespace, soupwise statistics')
% zlim([0 0.25])
% zlim([45 130])
% scatter3(log2(list1),log2(list2),log2(list3))
hold off
%%
ax=gca;
z=ax.ZLim(2);
% z=max(:)qwQ
zstep=2500;

for i=1:ceil(z/zstep)
    zlim([i-1 i]*zstep)
    disp(sprintf('%d~%d of %d',(i-1)*zstep,i*zstep,z))
    pause
    
end