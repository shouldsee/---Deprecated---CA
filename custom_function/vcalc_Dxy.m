function [ D ] = vcalc_Dxy( X,Y )
D=[];
elenum=max(X(:),max(Y(:)))+1;
for i=1:size(X,1);
   out=calc_Dxy(X(i,:),Y,elenum); 
   [m,n]=size(out);
   out=permute(out,[3 1 2]);
   D(i,1:m,1:n)=out;
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

end
end

