function [ cells] = ttorus( cells )
global x y
[a,b,c]=size(cells);
cells(a-1,:,:)=cells(1,:,:)+cells(a-1,:,:);
cells(2,:,:)=cells(end,:,:)+cells(2,:,:);
cells(:,b-1,:)=cells(:,1,:)+cells(:,b-1,:);
cells(:,2,:)=cells(:,end,:)+cells(:,2,:); 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

