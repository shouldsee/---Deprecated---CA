function [ B ] = diffcomb( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
elenum=max(A(:))+1;
A0=squeeze(A(1:end-1,:));
A1=squeeze(A(2:end,:));
B=A0+A1.*elenum;

end

