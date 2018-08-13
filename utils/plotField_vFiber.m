function plotField_vFiber(wf,wf0, k, info)

N = info.N;
xvals = info.xvals; 
yvals = info.yvals;    
lambdaOverD = info.lambdaOverD;

Gu = wf-wf0;

fiberDiam_pix = info.fiberDiam_pix;
x_fib_pix = info.x_fib_pix;
y_fib_pix = info.y_fib_pix;
wf_trim = wf(N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2,N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2);
wf_trim = wf_trim/sqrt(info.normI);
Gu_trim = Gu(N/2+y_fib_pix-fiberDiam_pix/2+1:N/2+y_fib_pix+fiberDiam_pix/2,N/2+x_fib_pix-fiberDiam_pix/2+1:N/2+x_fib_pix+fiberDiam_pix/2);
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