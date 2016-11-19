

rulelist=table2cell(T1(:,1));
% tp=table2array(T1(:,2:14));
% x=tp(:,1);
% y=tp(:,2);
% z=tp(:,3);
% ind=y-x<-0.3 & x<3 & y<1.6 & z<0.5;
% T1=T1(ind,:);
stats=table2array(T1(:,2:end-1));
estep=stats(:,end);

% [~,ind]=sort(stats(:,6),'descend');
% rulelist=rulelist(ind,:);
% stats=stats(ind,:);

% list1=stats(:,7);
% list2=stats(:,8);
% list3=stats(:,9);

% lists=reshape(stats(:,),[],3)
x=repmat(stats(:,1),2,1);
y=repmat(stats(:,2),2,1);
% z=repmat(stats(:,3),3,1);

% list1=reshape(stats(:,[4 7 10]),[],1)*0.2;
% list2=reshape(stats(:,[4 7 10]+1),[],1)*1;
% list3=reshape(stats(:,[4 7 10]+2),[],1)*500;
% id0=[4 7 10];
id0=[3 5];
list1=reshape(stats(:,id0),[],1);
list2=reshape(stats(:,id0+1),[],1);
% list3=reshape(stats(:,[4 7 10]+2),[],1);
% list1=list1./sqrt(sum([list1,list2,list3],2));
% list2=list2./sqrt(sum([list1,list2,list3],2));
list1=list1./sqrt(sum([list1,list2],2));
list2=list2./sqrt(sum([list1,list2],2));
% list3=list3./sqrt(sum([list1,list2,list3],2));
% list1=list1/4;
% list2=list2
rulelist=repmat(rulelist,3,1);
% list3=stats(:,9);

% list4=stats(:,6);
% ind=find(stats(:,8)<100 | list3<0.30|list2<0.1);
% x=x/2.5;
% y=y/2.5;
% z=z/2.5;
% ind=find(repmat(estep,2,1)<76|list1>3.5|(abs(list1)+abs(list2)<0.5));
ind=find(repmat(estep,2,1)<76 | x>3.5|y<5);
% ind=find(ind | (abs(list1)+abs(list2)<0.5));
x(ind)=-10;
% quiver3(x,y,z,list1,list2,list3,0,'MaxHeadSize',0.01);
quiver3(x,y,repmat(1:(length(x)/2),1,2)',list1,list2,ones(size(x))*0,0,'MaxHeadSize',0.01);

% iind=find(~ismember(1:length(rulelist),ind));
% filename='filtered';
% [rulename,uind]=unique(rulelist(iind));
% l2=list2(iind);
% l4=list4(iind);
% T=table(rulename,l2(uind),l4(uind));
% update_file
% writetable(table(rulename),'filtered.xls')
% ind=find(list2>1.0 & list3>100);
% rulename=unique(rulelist(ind));


% scatter3(list1,list2,list3,3,'x');
% gscatter(list2,list3,g,colormap(gca),'x',5);
