soname='so_SI';
load(soname,'so')
%%

% load visual_test
% load weirdo
% load ambig2
% load visual_ambig
load ambig_und
% load visual_chaos
% load life
% load visual_complex
figure(2)
clf

scatter(0,0,3,'x')
xlim([0 4])
ylim([0 2])
hold on
rulename=unique(rulename);
for i=1:length(rulename)
    rule=rulename{i}
try
rgroup=so.MapObj(rule);
end
id=find(so.group==rgroup);
% scatter(so.stats(id,1),so.stats(id,2),3,'x')
plot3(so.stats(id,1),so.stats(id,2),1:length(so.stats(id,1)))

hold on
drawnow
pause
end