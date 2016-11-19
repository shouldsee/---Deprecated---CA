
% clear
% load('B045678','rulename');

name='nerd';
load(name,'rulename');
% load('ambig_und','rulename');
% rulename=table2cell(readtable('nerd.xls'));

% excess=temp4
% excess=unique(excess);
% rulename=excess;
ruletemp=importrule(rulename);
restart
corner='A1';
% filename='B045678_test1';
filename=['flow_' name];
stats=stats(1:soupmax*length(rulename),:);
rulelist=rulelist(1:soupmax*length(rulename),:);
souplist=souplist(1:soupmax*length(rulename),:);
clear rule
save(name)
print2file
dirname='gallery/test/';
%%
% scanthrough
% save visual_complex