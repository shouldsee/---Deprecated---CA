function [ cells] = torus( cells )
global x y
cells(1,:)=cells(x(numel(x)),:);
cells(end,:)=cells(x(1),:);
cells(:,1)=cells(:,y(numel(y)));
cells(:,end)=cells(:,y(1)); 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

