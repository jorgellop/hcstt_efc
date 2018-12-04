% test_addAxisLabels4Paper
%
% Read in image and add axis labels
%
% Jorge Llop Sayson - Oct 18, 2018
close all

path = '/Users/jllopsay/Documents/GitHub/hcstt_efc/testing/';
namefile_arr = ["CamImage_flatDMRegEFC","CamImage_flatDMSMF","CamImage_flatDMBB",...
    "CamImage_finalRegEFC","CamImage_finalEFCSMF","CamImage_finalBroadband"];
for II=1:numel(namefile_arr)
    namefile = namefile_arr{II};
    im = imread([path,namefile,'.png']);
    fig0 = figure('visible','off','color','w');
    % fig0 = figure('visible','off');
    imagesc(im)
    axis image
    set(gca,'xtick',[],'ytick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',16)
%     set(gca,'visible','off')
    % hAx = axes;
    % set(gca, 'XColor','w','YColor','w')
    xlabel('{\itx (pixels)}')
    ylabel('{\ity (pixels)}')
    box off
    figure(fig0)
    export_fig([path,namefile,'_new.png'],'-r300');
end
