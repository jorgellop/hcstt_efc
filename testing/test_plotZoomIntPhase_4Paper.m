path_regEFC = '/Users/jllopsay/Documents/MATLAB/FiberCoupling_HCST/output/EFC_contrastComparison_Oct0418_regEFC_broadband_errormap2_12act/';
path_SMF = '/Users/jllopsay/Documents/MATLAB/FiberCoupling_HCST/output/EFC_contrastComparison_Oct0418_fiber_broadband_errormap2_12act/';
close all

lambdaOverD = 8;
fontsz = 16;
fiberDiam_pix = 1.4*8;
fig0 = figure('visible','off','color','w','pos',[10 10 1100 530]);
for II=1:4
    subplot(2,4,II)
    load([path_regEFC,'IntPhaseOverFibxel_it',num2str(II)])
    imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
        [y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
        ph_wf,[-pi pi]);
    axis image;
    set(gca,'FontSize',fontsz);
    titl = ['Iteration ',num2str(II)];
    title(titl,'FontSize', fontsz);
    colormap('hsv');
    set(gca,'tickdir','out','TickLength',[0.03, 0.01]) 
    xlabel('\itAngular separation (\lambda/D)','FontSize',fontsz);
    if II==1; ylabel({'\fontsize{21}Conventional EFC';'\fontsize{18}\itAngular separation (\lambda/D)'});end
end
hp4 = get(subplot(2,4,4),'Position')
% h = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.065  0.035  hp4(4)-0.132]);
% ylabel(h, 'Log(Raw Contrast)')
for II=1:4
    subplot(2,4,II+4)
    load([path_SMF,'IntPhaseOverFibxel_it',num2str(II)])
    imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
        [y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
        ph_wf,[-pi pi]);
    axis image;
    set(gca,'FontSize',fontsz)
    titl = ['Iteration ',num2str(II)];
    title(titl,'FontSize', fontsz);
    xlabel('\itAngular separation (\lambda/D)','FontSize',fontsz);
    set(gca,'tickdir','out','TickLength',[0.03, 0.01]) 

    colormap('hsv');
    if II==1; ylabel({'\fontsize{21}Fiber-based EFC';'\fontsize{18}\itAngular separation (\lambda/D)'});end
end
hp8 = get(subplot(2,4,8),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp8(2)  0.035  hp4(4)*2+0.132]);
ylabel(h, 'Phase (Rad)')
figure(fig0)
export_fig('ZoomIn_Ph_Oct4.png','-r300');

%%
fig0 = figure('visible','off','color','w','pos',[10 10 1100 530]);
for II=1:4
    subplot(2,4,II)
    load([path_regEFC,'IntPhaseOverFibxel_it',num2str(II)])
    imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
        [y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
        log10(wf_trim),[-7.5 -4.5]);
    axis image;
    set(gca,'FontSize',fontsz);
    titl = ['Iteration ',num2str(II)];
    title(titl,'FontSize', fontsz);
    xlabel('\itAngular separation (\lambda/D)','FontSize',fontsz);
    set(gca,'tickdir','out','TickLength',[0.03, 0.01]) 
    if II==1; ylabel({'\fontsize{21}Conventional EFC';'\fontsize{18}\itAngular separation (\lambda/D)'});end
end
hp4 = get(subplot(2,4,4),'Position')
% h = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.065  0.035  hp4(4)-0.132]);
% ylabel(h, 'Log(Raw Contrast)')
for II=1:4
    subplot(2,4,II+4)
    load([path_SMF,'IntPhaseOverFibxel_it',num2str(II)])
    imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
        [y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
        log10(wf_trim),[-7.5 -4.5]);
    axis image;
    set(gca,'FontSize',fontsz)
    titl = ['Iteration ',num2str(II)];
    title(titl,'FontSize', fontsz);
    xlabel('\itAngular separation (\lambda/D)','FontSize',fontsz);
    set(gca,'tickdir','out','TickLength',[0.03, 0.01]) 
    if II==1; ylabel({'\fontsize{21}Fiber-based EFC';'\fontsize{18}\itAngular separation (\lambda/D)'});end
end
hp8 = get(subplot(2,4,8),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp8(2)  0.035  hp4(4)*2+0.132]);
ylabel(h, 'Log(Raw Contrast)')
figure(fig0)
export_fig('ZoomIn_Int_Oct4.png','-r300');


