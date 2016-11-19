function [ E] = entropise( X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X=X/sum(X(:));
E=-X.*log2(X);
E(isnan(E))=0;
end

