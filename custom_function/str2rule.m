function[rule]=str2rule(str);

s=strsplit(str,'/');
Bnum=str2num(strip(s{1},'B')')';
Snum=str2num(strip(s{2},'S')')';
ruleB=[];
ruleS=[];
for i=0:8
    ruleB=[ruleB ismember(i,Bnum)];
    ruleS=[ruleS ismember(i,Snum)];
end
rule=[ruleB ruleS];
end