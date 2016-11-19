function [Y ] = findneighbor( X,Vset)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[m n]=size(Vset);

if ndims(X)==2;
    Y=repmat(X,1,1,m);
else
    Y=X;
end
% Y=zeros([size(X),m]);
for i=1:m
    Y(:,:,i)=circshift(Y(:,:,i),Vset(i,:));

end
end

