function plotIntPhase_vFiber_broadband(wf,wf0, k, info)

N = info.N;
xvals = info.xvals; 
yvals = info.yvals;    
lambdaOverD = info.lambdaOverD;
numOfWavelengths = info.numOfWavelengths;
lambda0 = info.lambda0;
lam_arr = info.lam_arr;
iGu_BB = zeros(N);
iWf_BB = zeros(N);

Gu = wf-wf0;

for II = 1:numOfWavelengths
    bm_lam = wf(:,:,II);
    Gu_lam = Gu(:,:,II);
    lam = lam_arr(II);
    lam_frac = lam/lambda0;
    FP = (lambda0/lam)*bm_lam;
    Gu_FP = (lambda0/lam)*Gu_lam;
    FP_lam = interp2(xvals,yvals,FP,...
             xvals/lam_frac,yvals/lam_frac,'linear',1e-30);
    Gu_FP_lam = interp2(xvals,yvals,Gu_FP,...
             xvals/lam_frac,yvals/lam_frac,'linear',1e-30);
    iWf_BB = iWf_BB + abs(FP_lam).^2/info.normI/numOfWavelengths;
    iGu_BB = iGu_BB + abs(Gu_FP_lam).^2/info.normI/numOfWavelengths;
end
fiberDiam_pix = info.fiberDiam_pix;
fiberDiam_pix = fiberDiam_pix(round(numOfWavelengths/2));
x_fib_pix = info.x_fib_pix;
y_fib_pix = info.y_fib_pix;
x_fib_pix = x_fib_pix(round(numOfWavelengths/2));
y_fib_pix = y_fib_pix(round(numOfWavelengths/2));
wf_trim = iWf_BB(N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2,N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2);
% wf_trim = wf_trim/sqrt(info.normI);
Gu_trim = iGu_BB(N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2,N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2);
% Gu_trim = Gu_trim/sqrt(info.normI);
%Calculate the phase
Gu0 = Gu(N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2,N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2, round(numOfWavelengths/2));
ph_Gu = angle(Gu0);
wffib = wf(N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2,N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2, round(numOfWavelengths/2));
ph_wf = angle(wffib);

full_namemat = [info.outDir,'IntPhaseOverFibxel_it',num2str(k),'.mat'];
full_path = [info.outDir,'IntPhaseOverFibxel_it',num2str(k),'.png'];

fig0 = figure('visible','off','color','w');
subplot(2,2,3)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
    [y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,ph_wf);
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Image: Phase over Fibxel, iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,1)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
    [y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,...
    log10(wf_trim),[-9 info.maxint]);
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Image: Intensity over Fibxel (log scale), iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,4)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,ph_Gu);
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Gu: Phase over Fibxel, iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,2)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,log10(Gu_trim),[info.minint info.maxint]);
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Gu: Intensity over Fibxel (log scale), iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);


export_fig(full_path,'-r300');
save(full_namemat,'wf_trim','ph_wf','x_fib_pix','y_fib_pix');
close(fig0);

end