function[rms,rulename]=single_mutate(rule);
if class(rule)=='char';
    r=importrule({rule});
    r=[r{1},r{2}];
elseif class(rule)=='cell';
    r=importrule(rule);
    r=[r{1},r{2}];
else
    r=rule;
end

rms=[];
rulename=[];
for i=0:numel(r);
    rm=r;
    if i==0;
    else
    rm(i)=1-rm(i);
    end
    rms=[rms;{rm(1:9),rm(10:18)}];
    rulename=[rulename;{rule2name({rm(1:9),rm(10:18)})}];
end
end