function plotFocusAmp_vFiberCoupling_broadband( bm, titl, full_path, info )
%plotFocusAmp Makes a nice plot of the pupil amplitude
%   Detailed explanation goes here
    N = info.N;
    xvals = info.xvals; 
    yvals = info.yvals;    
    lambdaOverD =  info.lambdaOverD;
    numOfWavelengths = info.numOfWavelengths;
    lambda0 = info.lambda0;
    lam_arr = info.lam_arr;
    iPSF_BB = zeros(N);
    for II = 1:numOfWavelengths
        bm_lam = bm(:,:,II);
        lam = lam_arr(II);
        lam_frac = lam/lambda0;
        FP = (lambda0/lam)*bm_lam;
        FP_lam = interp2(xvals,yvals,FP,...
                 xvals/lam_frac,yvals/lam_frac,'linear',1e-30);
        iPSF_BB = iPSF_BB + abs(FP_lam).^2/info.normI/numOfWavelengths;
    end
    ma = max(max(log10(iPSF_BB)));
    fig0 = figure('visible','off','color','w');
    imagesc(xvals,yvals,log10(iPSF_BB), [-10 ma]);
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

