function [ mx,my,m] = makeM( B )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[FIRX,FIRY]=ndgrid(-1:1,-1:1);
FIRM=ones(3,3);

m=conv2(B,FIRM,'same');
m=torus(m);
mx=conv2(B,FIRX,'same');
mx=torus(mx);
mx=mx./m;
% mx(isnan(mx))=0;

my=conv2(B,FIRY,'same');
my=torus(my);
my=my./m;
% my(isnan(my))=0;

end

