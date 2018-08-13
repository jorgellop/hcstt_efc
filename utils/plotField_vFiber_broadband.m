% Not using this function...

function plotField_vFiber_broadband(wf,wf0, k, info)

N = info.N;
xvals = info.xvals; 
yvals = info.yvals;    
lambdaOverD = info.lambdaOverD;
lambda0 = info.lambda0;
Gu = wf-wf0;

fiberDiam_pix = info.fiberDiam_pix;
x_fib_pix = info.x_fib_pix;
y_fib_pix = info.y_fib_pix;
lam_arr = info.lam_arr;
numOfWavelengths = info.numOfWavelengths;
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
    iWf_BB = iWf_BB + FP_lam/numOfWavelengths;
    iGu_BB = iGu_BB + Gu_FP_lam/numOfWavelengths;
end

wf_trim = iWf_BB(N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2,N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2);
wf_trim = wf_trim/sqrt(info.normI);
Gu_trim = iGu_BB(N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2,N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2);
Gu_trim = Gu_trim/sqrt(info.normI);

full_path = [info.outDir,'FieldOverFibxel_it',num2str(k),'.png'];

fig0 = figure('visible','off','color','w');
subplot(2,2,1)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,real(wf_trim),[info.minre info.maxre]);
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Image: Re(wf) over Fibxel (log scale), iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,2)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,real(Gu_trim));
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Gu: Re(Gu) over Fibxel (log scale), iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,3)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,imag(wf_trim),[info.minim info.maxim]);
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Image: Im(wf) over Fibxel (log scale), iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,4)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,imag(Gu_trim));
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Gu: Im(Gu) over Fibxel (log scale), iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

export_fig(full_path,'-r300');
close(fig0);

end