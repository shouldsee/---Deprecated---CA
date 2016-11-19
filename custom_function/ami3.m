%%% auto MI with respect to a 3d matrix
S_inputo=reshape(LS_inputo,size(LS_inputo,1),30,30);
X=S_inputo(1:10,:,:);
elenumtp=max(X(:))+1;
[a b c]=size(X);
as=1:2;
bs=1:2;
cs=1:2;
% bs=floor(b/3):floor(2*b/3);
% cs=floor(c/3):floor(2*c/3);

[X1,X2,X3]=ndgrid(as,[1:6 25:b],[1:6 25:c]);
vecs=[X1(:),X2(:),X3(:)];
MIs=zeros(size(X));
Dx=hist(X(:),(1:elenumtp)-1);
Dx=Dx/sum(Dx);
Hx=sum(-Dx.*log2(Dx),'omitnan');
elenumtp=length(Dx)+1;
repX=repmat(X,[1,3,3]);
as=1:a;
% bs=1:b;
% cs=1:c;
% bs=(1:b)+floor(b/2);
% cs=(1:c)+floor(c/2);
for vi=1:length(vecs)
    vec=vecs(vi,:);
    FIRtp=1;FIRtp(vec(1),vec(2),vec(3))=elenumtp;FIRtp(1,1,1)=1; 
    Sxy=convn(repX,FIRtp,'same');
    Sxy=Sxy(:,bs,cs);
    Dxy=hist(Sxy(:),(1:elenumtp^2)-1);
    Dxy=Dxy/sum(Dxy);
    Hxy=sum(-Dxy.*log2(Dxy),'omitnan');
    
    MIs(vec(1),vec(2),vec(3))=2*Hx-Hxy;
%     FIRtp
end

%%
MIs=circshift(MIs,[0 15 15]);
%%
MIs=permute(MIs,[2 3 1]);
%%
% LMIs=reshape(MIs,30,[]);
% imagesc(LMIs)

ind=find(MIs>prctile(MIs(:),97.75));
[xs,ys,zs] = ind2sub(size(MIs), ind);
v=[xs,ys,zs];
scatter3(xs, ys, zs,10*MIs(ind));
% plot3()