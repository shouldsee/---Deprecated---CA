function [ruleB ruleS] =BSrule( B,S )

B=str2num(B')';
S=str2num(S')';
ruleB=[];
ruleS=[];
for i=0:8
    ruleB=[ruleB ismember(i,B)];
    ruleS=[ruleS ismember(i,S)];
end
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here


end

