function[CORR]=spatial_corr(A);
% [a b]=size(A);
h=fspecial('average',3);
B=conv2(A,h,'valid');
A_=A(2:end-1,2:end-1);
CORR=corr(A_(:),B(:));
end