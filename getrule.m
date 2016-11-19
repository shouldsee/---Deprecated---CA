if randrule==0
            r=ruletemp(rind,:);
             rule=[r{1},r{2}];
elseif randrule==1;
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
elseif randrule==2;
                r=rulecurr*(2.^((18:-1:1)-1)')+1;
                rule=str2num(dec2base(r,2,18)')';
                ruletemp{rind,1}=rulecurr(1:9);
                ruletemp{rind,2}=rulecurr(10:18);

                rulename{rind}=rule2name(rulecurr);
end
        rulecurr=rule;