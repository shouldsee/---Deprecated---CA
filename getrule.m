            if ~randrule
            r=ruletemp(rind,:);
             rule=[r{1},r{2}];
            else
          ktemp=k;
              k=8;
            rrule=1;
         update_rule;
           ruletemp(rind,:)={ruleB,ruleS};
         rulename{rind}=rstr;
         k=ktemp;
          rule=rulecurr;
    %     rstr
    %     pause
        end
        rulecurr=rule;