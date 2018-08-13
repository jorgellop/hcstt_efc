% test_matlabIndecesMatrix
%
% Test to see how indeces fall on matrix when assigning values via double
% for loop
%
% Jorge Llop - Aug 10, 2018
numreg = 7;
numgain = 5;
regval_arr = logspace(-8,-1,numreg);
gain_arr = linspace(0.1,4,numgain);
curr_reg_arr = zeros(1,numreg*numgain) + nan;
countreg = 1;
for III=1:numgain
    for II=1:numreg
    %     regval = 0.1; % To be determined
        if II==6
            if III == 2
                curr_reg_arr(countreg) = -1;
                countreg
            end
        end
        countreg = countreg+1;
    end
end
[mi,ind_min] = min(curr_reg_arr);
[ind_minreg,ind_mingain] = ind2sub([numreg,numgain],ind_min)
% regval = regval_arr(ind_minreg);
% gainval = gain_arr(ind_mingain);
