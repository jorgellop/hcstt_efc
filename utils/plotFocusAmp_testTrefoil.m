function plotFocusAmp_testTrefoil( bm, titl, full_path, info )
%plotFocusAmp Makes a nice plot of the pupil amplitude
%   Detailed explanation goes here

    xvals = info.xvals; 
    yvals = info.yvals;    
    
    fig0 = figure('visible','off','color','w');
    imagesc(xvals*bm.dx*1e3,yvals*bm.dx*1e3,fftshift(abs(bm.wf)), [0.0005 0.0055]);
    colorbar; 
    axis image;
    axis([-10 10 -10 10]*bm.wl*bm.fr*1e3);
    title(titl);
    xlabel('x (mm)');
    ylabel('y (mm)');
    text(5*bm.wl*bm.fr*1e3,9*bm.wl*bm.fr*1e3,['f/',num2str(bm.fr,'%.1f')],'Color','w');
    currticks = get(gca,'XTick');
    set(gca,'XTick',currticks,'YTick',currticks);
    set(gca,'TickDir','out');set(gca,'YDir','normal');
    set(fig0,'units', 'inches', 'Position', [0 0 5 5])
    export_fig(full_path,'-r300');
    close(fig0);
    
end

