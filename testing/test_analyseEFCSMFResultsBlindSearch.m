% test_analyseEFCSMFResultsBlindSearch
%
%
%
% Jorge Llop - Aug 28, 2018

folder = '/Users/jllopsay/Documents/GitHub/hcstt_efc/output/EFC_wFiber_LabDemonstration_FibervSearchPos_35Positions_Aug28/';
filelist = dir([folder,'data*.mat']);
% filelist_arr = filelist(:).name;
numfile = numel(filelist);
for II=1:numfile
    [actxc_fib,ang_fib] = hcstt_test_getSpatialFreqOfFiber(II);
    actxc_fib_arr(II) = actxc_fib;
    ang_fib_arr(II) = ang_fib;
    load([folder,filelist(II).name])
    supp_arr = int_in_DH(1)./int_in_DH;
    maxsupp(II) = max(supp_arr);
    supp_it2(II) = int_in_DH(1)/int_in_DH(2);
    supp_lastit(II) = int_in_DH(1)/int_in_DH(numel(int_in_DH));
end
figure(301)
scatter(actxc_fib_arr,ang_fib_arr,maxsupp,maxsupp,'filled')
colorbar
title('Max Suppression')
figure(302)
scatter(actxc_fib_arr,ang_fib_arr,supp_it2,supp_it2,'filled')
colorbar
title('Suppression at iteration 2')
figure(303)
scatter(actxc_fib_arr,ang_fib_arr,supp_lastit,supp_lastit,'filled')
colorbar
title('Last Iteration Suppression')

