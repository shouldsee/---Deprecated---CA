function[Y]=circenum(X,vec)
% X=S_pop;
n=size(X)./vec;
n(isnan(n)|isinf(n))=0;
n=max(n);

Y=zeros([n,size(X)]);
for i=0:n-1
    Y(i+1,:)=reshape(circshift(X,padarray(i*vec,ndims(X)-length(vec),0,'post')),1,[]);   
end