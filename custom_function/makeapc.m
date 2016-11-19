function [ apc ] = makeapc( input,varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
update=varargin{1};
window=6;
        apc=[];
        l=length(input);
     for i=1:l-window
        wd=input(i:i+window);
        a=autocorr(input(i:i+window),window-1);
        apc=[apc;a];
       
if update==1;
       figure(11)

     waterfall(apc);
     title('apc waterfall');
      end
end

