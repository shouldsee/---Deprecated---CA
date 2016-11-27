function[rstr]=rule2name(rule)
ruleB=rule{1};
ruleS=rule{2};
B=find(ruleB==1)-1;
S=find(ruleS==1)-1;
rstr=sprintf('B%s/S%s',num2str(B')',num2str(S')');
end