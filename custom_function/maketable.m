function[values]=maketable(A,m)

% A=[1 6;36 216];
% m=6;
l=length(A(:));
%%
left=[];
right=[];
right1=[];
for i=1:l
    Var=['x' num2str(i)];
    left=[left Var ' '];
    right=[right '0:' num2str(m-1) ','];
    right1=[right1 Var '(:) ' ];
end
right(end)=[];
cd=sprintf('[%s]=ndgrid(%s);',left,right);
eval(cd)
cd=sprintf('X=[%s];',right1);
eval(cd)
size(X)

%%
enumer=zeros(2*l,l);
value=1;
values=zeros(2*l,1);
for i=1:length(X)
    vec=X(i,:);
    rvec=vec(end:-1:1);
    for j=(1:l)-1
        enumer(2*j+1,:)=circshift(vec,[0 j]);
        enumer(2*j+2,:)=circshift(rvec,[0 j]);
    end
    uenumer=unique(enumer,'rows');
    ind=find(ismember(X(1:i-1,:),uenumer,'rows'));
    if isempty(ind);
        values(i)=value;
        value=value+1;
    else
        values(i)=values(min(ind));
    end
        
end
%%
v0=sum(bsxfun(@times,X,A(:)'),2);
values=values(v0+1);