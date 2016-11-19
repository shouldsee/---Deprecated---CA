function[Y]=stack(X)
% X=LS_input_old;
%%
elenum.X=max(X(:))+1;
Y=zeros([size(X,1) size(X,1) size(X,2)]);
%%
m=size(X,1);
for i=1:m
    for j=i:m;
%         x=permute(X(i:j,:),[1,3,2]);
        x=X(i:j,:);
        etemp=max(x(:))+1;
        y=sum(bsxfun(@times,etemp.^((i:j)-i)',x),1);
%         Y(i,j,:)=permute(y,[1 3 2]);
        Y(i,j,:)=y;
    end
end
