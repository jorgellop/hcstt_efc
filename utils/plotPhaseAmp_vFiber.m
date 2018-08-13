function plotPhaseAmp_vFiber(wf,wf0, k, info)

N = info.N;
xvals = info.xvals; 
yvals = info.yvals;    
lambdaOverD = info.lambdaOverD;

Gu = wf-wf0;

fiberDiam_pix = info.fiberDiam_pix;
x_fib_pix = info.x_fib_pix;
y_fib_pix = info.y_fib_pix;
wf_trim = wf(N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2,N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2);
wf_trim = wf_trim/sqrt(info.normI);
Gu_trim = Gu(N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2,N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2);
Gu_trim = Gu_trim/sqrt(info.normI);

full_path = [info.outDir,'IntPhaseOverFibxel_it',num2str(k),'.png'];

fig0 = figure('visible','off','color','w');
subplot(2,2,3)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,angle(wf_trim));
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Image: Phase over Fibxel, iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,1)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,log10(abs(wf_trim).^2));
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Image: Intensity over Fibxel (log scale), iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,4)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,angle(Gu_trim));
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Gu: Phase over Fibxel, iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

subplot(2,2,2)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,log10(abs(Gu_trim).^2));
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Gu: Intensity over Fibxel (log scale), iteration ',num2str(k)];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);


export_fig(full_path,'-r300');
close(fig0);

end