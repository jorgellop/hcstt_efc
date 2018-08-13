function hcstt_test_plotModelImage(im, full_path, info )
%plotFocusAmp Makes a nice plot of the pupil amplitude
%   Detailed explanation goes here

%     xvals = info.xvals; 
%     yvals = info.yvals;    
    N = info.N;
    normI = info.normPower;
    fig0 = figure('visible','off','color','w');
%     imagesc(im,[-0.005 0.005]);
%     imn = im/max(im(:));
    imn = im;
%     imagesc(log10((abs(imn(N/2-200:N/2+200,N/2-200:N/2+200)).^2)), [-10 0])
    imagesc(log10((abs(imn(N/2-200:N/2+200,N/2-200:N/2+200)).^2)))
    colorbar; 
    axis image;
%     axis([-10 10 -10 10]*bm.wl*bm.fr*1e3);
    title('HCSTT Camera Image');
    xlabel('x (pix)');
    ylabel('y (pix)');
%     text(5*bm.wl*bm.fr*1e3,9*bm.wl*bm.fr*1e3,['f/',num2str(bm.fr,'%.1f')],'Color','w');
    currticks = get(gca,'XTick');
    set(gca,'XTick',currticks,'YTick',currticks);
    set(gca,'TickDir','out');set(gca,'YDir','normal');
    set(fig0,'units', 'inches', 'Position', [0 0 5 5])
    export_fig(full_path,'-r300');
    figure(1)
%     close(fig0);
    
end

