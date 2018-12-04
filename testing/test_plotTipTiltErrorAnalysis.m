% test_plotTipTiltErrorAnalysis
%
% Read data computed previously and plot a 4*2 sublot plot
%
% Jorge Llop - Oct 5, 2015
clear all;close all;
fontsz = 15;
path = '/Users/jllopsay/Documents/MATLAB/HCSTR_fiber_proper_model/output/HCST_FIU_broadband12act_loadAndRuntest_Oct04_EFCwDM1_VC6_N1024_circ_LSout0.9_lam6.5e-07_BW10_surf_error_PSDfit/';
load([path,'tiptiltErrorAnalysis'])
load([path,'tiptiltErrorAnalysis_direction2And4'])

direction_arr = [0,90,180,270,45,90+45,180+45,270+45];
numtry = 5;
maxOff = 0.02;
vect = linspace(0,maxOff,numtry)/8*100;
% vect = linspace(0,1,numtry);
szfig = [1200 530];
fig0 = figure('visible','off','color','w','pos',[10 10 szfig(1) szfig(2)]);
numOfWavelengths=31;
for II=1:4
    for k=1:5
        handle=subtightplot(2,4,II);
        if II==2; plot(lam_arr*1e9,log10(S_in_DH_fib_mat2_direction2And4(k,:,1)),'LineWidth',3);
        elseif II==4; plot(lam_arr*1e9,log10(S_in_DH_fib_mat2_direction2And4(k,:,2)),'LineWidth',3);
        else
            plot(lam_arr*1e9,log10(S_in_DH_fib_mat2(k,:,II)),'LineWidth',3);
        end
        hold on
        
        
%         hold on
%         axes('pos',[10 10 1 10])
%         imshow(['arrowOffset',num2str(k),'.png'])
        % p(1).LineWidth = 3;
        % p(2).LineWidth = 3;
        % p(1).Color = 'r';
        % p(2).Color = 'b';
        %     hold on
        %     plot(1:k,log10(coupl_MMF_in_DH))
%         xlabel('Wavelength [nm]','fontsize',fontsz)
        set(gca,'XTickLabel',[])
        if II==1; ylabel('log(SMF Raw Contrast)','fontsize',fontsz); else
            set(gca,'YTickLabel',[]); end
        set(gca,'FontSize',fontsz);
        axis([ lam_arr(1)*1e9 lam_arr(numOfWavelengths)*1e9 -12 -5])
        h = vline(650,'r','\lambda_0');
        h = vline(617.5,'b','');
        h = vline(682.5,'b','');

%         ha =gca;
%         haPos = get(handle,'position');
%         ha2=axes('position',[0.7, 0.5, .02,.02,]);
%         test=imread(['arrowOffset',num2str(II),'.png']);
%         image(te st)
%         set(ha2,'handlevisibility','off','visible','off')  
    end
    hold off
    lgd = legend('No displacement',[num2str(vect(2)),'%\lambda/D'],[num2str(vect(3)),'%\lambda/D'],[num2str(vect(4)),'%\lambda/D']);%,'Location','northeastoutside');
    title(lgd,{['Displacement direction: '],[num2str(direction_arr(II)),'^{\circ}  ']})
%     title(['Displacement direction: ',num2str(direction_arr(II)),'^{\circ}'],'fontsize',16)
    lgd.FontSize = 12;
    lgdpos = lgd.Position;
    lgdpos(1) = lgdpos(1)+0.008;
    lgdpos(2) = lgdpos(2)+0.0181;
    lgd.Position = lgdpos;
    haPos = get(handle,'position');
    impos = lgdpos;
    impos(1) = impos(1)-0.0937; 
%     if II==1; impos(1) = impos(1)-0.0605; end
%     if II==2; impos(1) = impos(1)-0.0565; end
%     if II==3 || II==4; impos(1) = impos(1)-0.0525; end
    impos(2) = impos(2);
    impos(3) = impos(4)*szfig(2)/szfig(1);
    ha2=axes('position',impos);
    test=imread(['arrowOffset',num2str(II),'.png']);
    image(test)
    set(ha2,'handlevisibility','off','visible','off')  
    yticks([-11:2:-7])
figure(fig0)
end
for II=5:8
    for k=1:5
        handle=subtightplot(2,4,II);
%         handle=subplot(2,4,II);
        xlabel('Wavelength [nm]','fontsize',fontsz)
        if II==2; plot(lam_arr*1e9,log10(S_in_DH_fib_mat2_direction2And4(k,:,1)),'LineWidth',3);
        elseif II==4; plot(lam_arr*1e9,log10(S_in_DH_fib_mat2_direction2And4(k,:,2)),'LineWidth',3);
        else
            plot(lam_arr*1e9,log10(S_in_DH_fib_mat2(k,:,II)),'LineWidth',3);
        end
        hold on
        
        
%         hold on
%         axes('pos',[10 10 1 10])
%         imshow(['arrowOffset',num2str(k),'.png'])
        % p(1).LineWidth = 3;
        % p(2).LineWidth = 3;
        % p(1).Color = 'r';
        % p(2).Color = 'b';
        %     hold on
        %     plot(1:k,log10(coupl_MMF_in_DH))
        
        if II==5; ylabel('log(SMF Raw Contrast)','fontsize',fontsz); else
            set(gca,'YTickLabel',[]); end
        set(gca,'FontSize',fontsz);
        axis([ lam_arr(1)*1e9 lam_arr(numOfWavelengths)*1e9 -12 -6])
        % lgd = legend('Coupling - EFC w/SMF');
%         lgd.FontSize = fontsz;
        h = vline(650,'r','\lambda_0');
        h = vline(617.5,'b','');
        h = vline(682.5,'b','');
%         ha =gca;
%         haPos = get(handle,'position');
%         ha2=axes('position',[0.7, 0.5, .02,.02,]);
%         test=imread(['arrowOffset',num2str(II),'.png']);
%         image(test)
%         set(ha2,'handlevisibility','off','visible','off')  
    end
    
    figure(fig0)
    lgd = legend('No displacement',[num2str(vect(2)),'%\lambda/D'],[num2str(vect(3)),'%\lambda/D'],[num2str(vect(4)),'%\lambda/D']);%,'Location','northeastoutside');
    title(lgd,{['Displacement direction: '],[num2str(direction_arr(II)),'^{\circ}  ']})
%     title(['Displacement direction: ',num2str(direction_arr(II)),'^{\circ}'],'fontsize',16)
    lgd.FontSize = 12;
    lgdpos = lgd.Position;
    lgdpos(1) = lgdpos(1)+0.008;
    lgdpos(2) = lgdpos(2)+0.0181;
    lgd.Position = lgdpos;
    haPos = get(handle,'position');
    impos = lgdpos;
    impos(1) = impos(1)-0.0937; 

%     if II==5; impos(1) = impos(1)-0.0565; end
%     if II==6; impos(1) = impos(1)-0.0525; end
%     if II==7 || II==8; impos(1) = impos(1)-0.0525; end
    impos(2) = impos(2);
    impos(3) = impos(4)*szfig(2)/szfig(1);
    ha2=axes('position',impos);
    test=imread(['arrowOffset',num2str(II),'.png']);
    image(test)
    set(ha2,'handlevisibility','off','visible','off')  
    yticks([-11:2:-7])
figure(fig0)

end
hold off
export_fig('JitterAnaylsis4Paperv3.png','-r300');
