function[d]=fdist(A,B)
%%A is a 1 by N vector
%%B is a M by N matrix

N=1:size(A,2);
M=size(B,1);
d=zeros(M,1);
for i=1:M
    d(i)=frechet(N',A',N',B(i,:)');
end
end
