% find a subset of configs specified by cond
% calculate uncertainty about outcome of cell
%%


%%
% clf

% out0=squeeze(rulecurr(convn(ndict,FIR.S_input,'valid')+1));
% HD0=entropise(hist(out0,0:1));
% H0=sum(HD0);

% rec2=[];
ind_all=(1:size(ndict,3));
inputs=squeeze(convn(ndict,FIR.S_input,'valid'));
lndicts=char();
for i=1:512
    lndicts(i,:)=dec2base(i-1,2,9);
end
lndict=reshape(ndict,9,512);
plndict=bsxfun(@times,lndict,2.^(0:9-1)');
plndict=plndict([4 1 2 3 5 6 7 8 9],:);
oplndict=plndict;
oplndict(1,:)=[];
for rind=rind0:size(ruletemp,1)
%%
getrule
% ind=[];
% 
% configset=ndict(:,:,ind);
outs=rulecurr(inputs+1);
% for
for maxng=1:17
    if maxng>=9
        cond=ndict(2,2,:)==1;
        categ=sum(plndict(1:maxng-8,:),1);

    else
        cond=(sum(sum(ndict,2),1)-ndict(2,2,:))<=maxng;
        categ=sum(oplndict(1:maxng,:),1);

    end
% ind_pos=find(cond);
% ind_neg=ind_all(~ismember(ind_all,ind_pos));
% linind(ind_pos)=1;
% linind(ind_neg)=0;
MI=calc_MI(categ,outs);
Hx=calc_Hx(categ,outs);
Hy=calc_Hy(categ,outs);
Hxy=calc_Hxy(categ,outs);
rec(maxng)=Hy;
end

% MI=calc_Hx(linind,outs);

% dout=hist(outset,0:1);
% HD=entropise(dout);

rec2(rind,:)=(rec);
% H2=sum(HD(:));
fprintf('rind=%d, rule %s, H=%.3d \n',rind,rulename{rind},max(rec));
subplot(2,2,1)
hold on 
plot(rec)
hold off
end
%%
subplot(2,2,1)
imagesc(rec2' )
subplot(2,2,2)
% plot(mean(rec2,2,'omitnan'))
plot(rind0:rind-1,rec2(rind0:rind-1,10))