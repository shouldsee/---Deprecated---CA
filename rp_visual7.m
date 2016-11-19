pre_init
preallocate
% ruletemp=importrule(rulename);
kmax=length(rulename)*soupmax;
if randrule
    kmax=100000000;
end

for k=k0:kstep:kmax

    initialise
    
if sum(S_popo(:))==0
        for i=1:length(lists)
           list=lists{i};
           eval(sprintf('%ss=[%ss;zeros(1,size(%ss,2))];',list,list,list))
        end
    
        continue
end

rind=floor((k-1)/soupmax)+1;    
    
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
    rind=floor((k-1)/soupmax)+1;
    suptitle(['mean-field p map' ' rulenumber='  num2str(rind) ' ' rulename{rind}])

    mfappx
    if store
        makeprint
    end
% rind=ceil((k)/soupmax)    
% if rind>size(rulename,1);
%     rind=size(rulename,1)
%     rind=1;
% end
if work==1
stats=[H_inputo FH_inputo];
so.stats=[so.stats;stats];
so.rulelist=[so.rulelist;repmat(rulename(rind,:),size(stats,1),1)];
so.label=[so.label;label*ones(size(stats,1),1)];
end
% subplot(2,2,2)
% hold off
% 
% subplot(2,2,3)
% hold off

% if store==1
%     dir='gallery/ents/';
%     s=rulename(rind);
%     s=strrep(s,'/','');
%     figname=sprintf('%s_soup%d.jpg',s{1},k);
% 
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
%     saveas(gcf,[dir figname])
% end
% inds=[inds;k];
% 
% Dv1o=detrend(v1o);
% timp_v1o=1-var(Dv1o)/var(v1o);
% Mv1o=mean(v1o);
% v1s=[v1s;Mv1o timp_v1o];
% fprintf('k=%d  Mv1o=%.3f timp_v1o=%.4f',k,Mv1o,timp_v1o)
% 
% Dv2o=detrend(v2o);
% timp_v2o=1-var(Dv2o')/var(v2o');
% Mv2o=mean(v2o);
% v2s=[v2s;Mv2o timp_v2o];


% [pc,score,latent,tsquare] = pca(D_inputo);

if update2==1 && stepnumber>20
makeplot
end
if pause2==1
pause
end
end