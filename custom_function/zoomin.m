function [ ] = zoomin( rect )
global px py n

px=repmat(linspace(px(rect(1,1),1),px(rect(2,1),1),n)',1,n);
py=repmat(linspace(py(1,rect(1,2)),py(1,rect(2,2)),n),n,1);

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

