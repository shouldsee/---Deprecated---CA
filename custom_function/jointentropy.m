function[JH]=jointentropy(A,B)
elenum=max(max(A(:)),max(B(:)))+1;


A=A.*elenum;
C=bsxfun(@plus,A,B);

if size(B,1)==1
D_C=hist(C,0:elenum^2-1);

else

% D_A=hist(A,0:elenum-1);
% D_B=hist(B',0:elenum-1)';
% D=bsxfun(@times,)
D_C=hist(C',0:elenum^2-1)';
end
D_C=bsxfun(@rdivide,D_C,sum(D_C,2));
% D_A=D_A./sum(D_A(:));
% D_B=D_B./sum(D_B(:));

% id=bsxfun(@times,D_A,permute(D_B,[2,3,1]));
% id=reshape(permute(id,[3,1,2]),[],elenum^2);
% MI=sum(-log2(id./D_C).*D_C,2,'omitnan');
JH=sum(-log2(D_C).*D_C,2,'omitnan');
end