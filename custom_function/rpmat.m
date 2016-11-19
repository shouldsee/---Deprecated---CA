function [ rp] = rpmat( A,thres,kernel)
% FIR=[1 0 -1; 0 0 0; -1 0 1];
dm=conv2(A,kernel,'same');
dm=abs(dm(2:end-1,2:end-1));
rp=dm>thres;

end

