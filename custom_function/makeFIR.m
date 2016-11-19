% S_inputo=reshape(LS_inputo,size(LS_inputo,1),30,30);
function[FIRtp]=makeFIR(elenum,Y)
% X=S_inputo;
% Y=zeros(2,2,2);
% l=2
% elenum=1;

%%
% pos=1:8;
yn=numel(Y);
ys=size(Y);
yd=ndims(Y);
pos=1:yn;
pvec=combvec(pos,pos);
pvec=unique(sort(pvec)','rows')';
% unique(pvec(2,:)-pvec(1,:))
% pvec=pvec(:,pvec(1,:)<pvec(2,:));
%%
FIRtp={};
for i=1:size(pvec,2);
    x=zeros(1,yn);
    
    for j=1:size(pvec,1);
        x(pvec(j,i))=elenum.^(j-1);
    end
    if numel(nonzeros(x))==2
    FIRtp=[FIRtp,reshape(x,ys)];
    end
end

end