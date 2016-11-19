function [ cells] = torus( cells )
global x y
[a,b]=size(cells);
cells(1,:)=cells(a-1,:);
cells(end,:)=cells(2,:);
cells(:,1)=cells(:,b-1);
cells(:,end)=cells(:,2); 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

