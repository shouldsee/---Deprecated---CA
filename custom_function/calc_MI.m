function[MI]=calc_MI(A,B)
% elenum=max(max(A(:)),max(B(:)))+1;
Hx=calc_Hx(A,B);
Hy=calc_Hy(A,B);
Hxy=calc_Hxy(A,B);
MI=Hx+Hy-Hxy;
end