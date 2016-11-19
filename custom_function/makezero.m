function [X ] = makezero( X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

X(isnan(X)|isinf(X))=0;
end

