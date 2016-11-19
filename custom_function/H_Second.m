function[d]=H_Second(A,B)

elenum=18;
% 
% A=A;
% B=LS_inputo;

Tinput=repmat(A,size(B,1),1)+B.*elenum;
if size(Tinput,1)==1
    D_Tinput=hist(Tinput,(1:elenum^2)-1);
else
D_Tinput=hist(Tinput',(1:elenum^2)-1)';
end
D_Tinput=D_Tinput./repmat(sum(D_Tinput,2),1,size(D_Tinput,2));
HD_Tinput=-D_Tinput.*log2(D_Tinput);
d=sum(HD_Tinput,2,'omitnan');

end

