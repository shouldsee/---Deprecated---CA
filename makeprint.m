%     dir='gallery/NHT1/';
    s=rulename(rind);
    s=strrep(s,'/','');
    figname=sprintf('soup%d_%s.jpg',k,s{1});
    get(0,'screenpixelsperinch');
    
    fig=gcf;
    fig.PaperPositionMode = 'auto';
%     fig.Position=[0 0 800 800];

    %     fig.PaperUnits='points';
    %     fig.InvertHardcopy = 'off';
    saveas(gcf,[galleryname batchname figname])