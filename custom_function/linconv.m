function[TT]=linconv(H,X)
% H=ones(3,3);
% H(2,2)=0;
% X=ones(nx,nx);
[mx nx]=size(X);
[hm hn]=size(H);
hmh=floor(hm/2);hnh=floor(hn/2);
T = convmtx2(H,mx,nx); 
os=[hm hn]+[mx nx]-1;
TT=reshape(full(T),[os mx*nx]);
TT=ttorus(TT);
st=size(TT);
TT=TT(2:st(1)-hmh,2:st(2)-hnh,:);
TT=reshape(TT,mx*nx,mx*nx);
end