function plotBBpsf4Paper( bm, full_path, info )
%plotBBpsf Makes a nice plot of the broadband psf
%   Detailed explanation goes here

    cmapmin = -9;
    cmapmax = -2;
    axismax = 10;
    axistickstep = 5;

    xvals = info.xvals; 
    yvals = info.yvals;    
    
    lambda0 = info.lambda0; 
    numOfWavelengths = info.numOfWavelengths;
    normI=info.normI;
    
    N = info.N;
    
%     fig0 = figure('visible','off','color','w');
%     imagesc(xvals*bm.dx*1e3,yvals*bm.dx*1e3,fftshift(abs(bm.wf)));
%     colorbar; 
%     axis image;
%     axis([-10 10 -10 10]*bm.wl*bm.fr*1e3);
%     title(titl);
%     xlabel('x (mm)');
%     ylabel('y (mm)');
%     text(5*bm.wl*bm.fr*1e3,9*bm.wl*bm.fr*1e3,['f/',num2str(bm.fr,'%.1f')],'Color','w');
%     currticks = get(gca,'XTick');
%     set(gca,'XTick',currticks,'YTick',currticks);
%     set(gca,'TickDir','out');set(gca,'YDir','normal');
%     set(fig0,'units', 'inches', 'Position', [0 0 5 5])
%     export_fig(full_path,'-r300');
%     close(fig0);

    iPSF_BB = zeros(N);

    for index = 1:numel(bm);

        bm_lam = bm(index);
        lam_frac = bm_lam.wl/lambda0;

        FP = (lambda0/bm_lam.wl)*fftshift(gather(bm_lam.wf));

        FP_lam = interp2(xvals,yvals,FP,...
                 xvals/lam_frac,yvals/lam_frac,'linear',1e-30);

        iPSF_BB = iPSF_BB + abs(FP_lam).^2/normI/numOfWavelengths;

    end

    bm_lam = bm(round(numOfWavelengths/2));
    
    fig0 = figure('visible','off','color','w');
    imagesc(xvals*bm_lam.dx/(bm_lam.wl*bm_lam.fr),yvals*bm_lam.dx/(bm_lam.wl*bm_lam.fr),log10(iPSF_BB));
    hcb = colorbar; %colormap(hot(256));
    ylabel(hcb, 'Log(Normalized Intensity)')
    axis image;
    caxis([cmapmin cmapmax]);
    axis([-axismax axismax -axismax axismax]);
%     title(titl);
    hx = xlabel('Angular coordinate (\lambda_0/{\itD})','FontSize', 10);
    hy = ylabel('Angular coordinate (\lambda_0/{\itD})','FontSize', 10);
    set(gca,'XTick',-axismax:axistickstep:axismax,'YTick',-axismax:axistickstep:axismax)
    set(gca,'TickDir','out');set(gca,'YDir','normal');
    set(gca,'FontSize', 16,...
        'TickDir','out',...
        'TickLength',[.02 .02]);
    set(fig0,'units', 'inches', 'Position', [0 0 5 5])
    export_fig(full_path,'-r300');
    close(fig0);

    outsize = 256;
    crop = N/2-outsize/2+1:N/2+outsize/2;
    fitswrite(iPSF_BB(crop,crop),[full_path,'.fits']);

    
end

