function plotFocusAmp_vFiberCoupling( bm, titl, full_path, info )
%plotFocusAmp Makes a nice plot of the pupil amplitude
%   Detailed explanation goes here

    xvals = info.xvals; 
    yvals = info.yvals;    
    lambdaOverD =  info.lambdaOverD;
    ma = max(max(log10((abs(bm).^2)/info.normI)));
    fig0 = figure('visible','off','color','w');
    imagesc(xvals,yvals,log10((abs(bm).^2)/info.normI), [-10 ma]);
    colorbar; 
    axis image;
    axis([-10 10 -10 10]*lambdaOverD);
    title(titl);
    xlabel('x (mm)');
    ylabel('y (mm)');
%     text(5*bm.wl*bm.fr*1e3,9*bm.wl*bm.fr*1e3,['f/',num2str(bm.fr,'%.1f')],'Color','w');
    currticks = get(gca,'XTick');
    set(gca,'XTick',currticks,'YTick',currticks);
    set(gca,'TickDir','out');set(gca,'YDir','normal');
    set(fig0,'units', 'inches', 'Position', [0 0 5 5])
    export_fig(full_path,'-r300');
    close(fig0);
    
end

