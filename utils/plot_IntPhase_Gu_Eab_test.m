function plot_IntPhase_Gu_Eab_test(wf,wf0, Eab_re, Eab_im, info)

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

Eab_max = info.Eab_max;
Eab_min = info.Eab_min;

full_path = [info.outDir,'IntPhase',filesep,'IntPhaseOverFibxel_Eabre',num2str(Eab_re),'_Eabim',num2str(Eab_im),'.png'];

fig0 = figure('visible','off','color','w');
pos1 = [0.4 0.6 0.225 0.225];
subplot('Position',pos1)
scatter(Eab_re,Eab_im,120,'filled','square')
grid on
axis equal
axis([Eab_min Eab_max Eab_min Eab_max])
set(gca,'FontSize',7)
titl = ['Sum(Eab x SMF mode), Field'];
title(titl,'FontSize', 7);
xlabel('Re{Sum(Eab x SMF mode)}','FontSize', 7);
ylabel('Im{Sum(Eab x SMF mode)}','FontSize', 7);

pos1 = [0.1 0.1 0.415 0.415];
subplot('Position',pos1)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,log10(abs(Gu_trim).^2),[-11.5 -6]);
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Gu: Int(Gu) over Fibxel (log scale'];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
ylabel('y (lam/D)','FontSize', 7);

pos1 = [0.5 0.1 0.415 0.415];
subplot('Position',pos1)
imagesc([x_fib_pix-fiberDiam_pix/2+1:x_fib_pix+fiberDiam_pix/2]/lambdaOverD,[y_fib_pix-fiberDiam_pix/2+1:y_fib_pix+fiberDiam_pix/2]/lambdaOverD,angle(Gu_trim),[-3 3]);
colorbar; 
axis image;
set(gca,'FontSize',7)
titl = ['Gu: Phase(Gu) over Fibxel'];
title(titl,'FontSize', 7);
xlabel('x (lam/D)','FontSize', 7);
set(gca,'YTick',[])
% ylabel('y (lam/D)','FontSize', 7);

export_fig(full_path,'-r300');
close(fig0);

end