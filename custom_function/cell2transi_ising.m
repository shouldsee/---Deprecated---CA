function [ transi ] =cell2transi_ising( cells,J)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sumnh=conv2((cells-0.5)*2,[0 1 0; 1 0 1; 0 1 0],'same');
transi=exp(J.*sumnh.*(0.5-cells)*2);

end

