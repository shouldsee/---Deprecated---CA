k=8
if ~exist('rule')
[a,b]=size(rulelist);
rule=[];
for i=1:a
    ruleB=rulelist{i,1};
    ruleS=rulelist{i,2};
    
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
    B=B(ismember(B,0:8));
    S=S(ismember(S,0:8));

    r=sprintf('B%s/S%s',num2str(B')',num2str(S')');
    rule=[rule;{r}];
end
end
%%
% logfile=[filename '.csv'];
% col={'rule','welist' ,'snrlist', 'mmlist'};
% col={'rule','welist' ,'snrlist', 'mmlist','popslist'};

% C=[];
% for i=col
% C=[C eval(i{1})];
% end
% T=table(rule,list1,list2,list3,list4);
if ~exist('header')
    header=true;
end

% if ~exist([filename '.xls'])
%     fclose(fopen([filename '.xls'], 'w'));
% end
T=table(rule,stats,souplist);
writetable(T,[filename '.xls'],'Range',corner, 'WriteVariableNames',header)
clear rule
% %%
% fileID=fopen(logfile,'w');
% b=length(col);
% 
% for i=1:b
% fprintf(fileID,'%s,',col{i});
% end
% fprintf(fileID,'\n',''); 
% spec={'%s ,';'%1.5f ,';'%1.5f ,';'%.5f,';'%.5f \n'};
% for i=1:a
%     for j=1:b
%         ent=T{i,j};
%         if j==1;
%         ent=ent{1};
%         end
%         fprintf(fileID,spec{j},ent);
% 
%     end
% end
        