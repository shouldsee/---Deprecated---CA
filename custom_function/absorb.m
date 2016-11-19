function [ cells] = absorb( cells )
global x y
[~,ma]=max(hist(cells(2,:),0:1));
cells(1,:)=ma-1;
[~,ma]=max(hist(cells(end-1,:),0:1));
cells(end,:)=ma-1;
[~,ma]=max(hist(cells(:,2),0:1));
cells(:,1)=ma-1;
[~,ma]=max(hist(cells(:,end-1),0:1));
cells(:,end)=ma-1; 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

