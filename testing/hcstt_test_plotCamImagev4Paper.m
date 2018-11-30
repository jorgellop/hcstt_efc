function hcstt_test_plotCamImagev4Paper(im,titlestr, full_path , sz)
%plotFocusAmp Makes a nice plot of the pupil amplitude
%   Detailed explanation goes here
% close;

    nx = sz(1);
    ny = sz(2);
    fig0 = figure('visible','off','color','w');
%     imagesc(im,[-0.005 0.005]);
%     imn = im/max(im(:));
    imagesc(im)
    h=colorbar; 
    ylabel(h, 'Intensity ','FontSize',12)
    axis image;
%     axis([-10 10 -10 10]*bm.wl*bm.fr*1e3);
    title(titlestr);
%     xlabel('x (pix)');
%     ylabel('y (pix)');
%     text(5*bm.wl*bm.fr*1e3,9*bm.wl*bm.fr*1e3,['f/',num2str(bm.fr,'%.1f')],'Color','w');
    currticks = get(gca,'XTick');
    set(gca,'XTick',currticks,'YTick',currticks);
    set(gca,'TickDir','out');set(gca,'YDir','normal');
    set(fig0,'units', 'inches', 'Position', [0 0 5 5])
    export_fig(full_path,'-r300');
%     close(fig0);
    
end

