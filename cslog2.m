function[X]=cslog2(X)

id=X~=0;
X(id)=log2(X(id));
X(~id)=0;
end