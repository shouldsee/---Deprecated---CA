    if rrule==1;    
    ruleB=str2num(dec2base(randi([0 2^(k+1)],1),2,k+1).')';
    ruleS=str2num(dec2base(randi([0 2^(k+1)],1),2,k+1).')';
    elseif import_rule==1;
    [ruleB,ruleS]=ruletemp{rulenumber,:};
    else
%     i=id(rulenumber);
    i=rulenumber;
    [ruleB,ruleS]=rulelist{i,:};
    end
    
    
    if k==9;
    B=find(ruleB==1)-1;
    S=find(ruleS==1)-2;
    elseif k==8
    B=find(ruleB==1)-1;
    S=find(ruleS==1)-1;
    else
    B=ruleB;
    S=ruleS;
    end
    disp(['ruleB=' num2str(B)]);
    disp(['ruleS=' num2str(S)]);
%     disp(['rulenum=' num2str(rulenumber)])
    soupnumber=0;
    rulecurr=[ruleB ruleS];
    rstr=sprintf('B%s/S%s',num2str(B')',num2str(S')');
