function [ dCS_input,base ] = autoMI( S_input )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
input=S_input;
CS_inputX=circenum(input,[1 0]);
CS_inputY=circenum(input,[0 1]);
CS_inputX=reshape(CS_inputX,size(CS_inputX,1),[]);
CS_inputY=reshape(CS_inputY,size(CS_inputY,1),[]);
dCS_input=pdist2(CS_inputX,CS_inputY,@calc_MI);
dCS_input=circshift(full(dCS_input),[15 15]);
base=dCS_input(15,15);
end

