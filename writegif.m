function [] = writegif( figname, gifname )
global gifcount figs
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

im=frame2im(getframe(figs.(figname).Parent.Parent));
[imind,cm]=rgb2ind(im,256);
gifname=strrep(gifname,'/','_');
if gifcount==1;
    imwrite(imind,cm,[gifname '.gif'],'gif','Loopcount',inf,'DelayTime',0.05);
    gifcount=0;
else
    imwrite(imind,cm,[gifname '.gif'],'gif','WriteMode','append','DelayTime',0.05);
end

end

