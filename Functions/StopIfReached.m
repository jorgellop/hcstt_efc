%function to tell Multistart to stop running if it reaches desired contrast
%early
%Yinzi Xin

function stop = StopIfReached(optimValues, state)
persistent foundLocal
stop = false;
switch state
    case 'init'
        foundLocal = []; % initialized as empty
    case 'iter'
        newf = optimValues.localsolution.Fval;
        eflag = optimValues.localsolution.Exitflag;
        % Now check if the exit flag is positive and
        % the new value differs from all others by at least 10^-10
        % If so, add the new value to the newf list
        if eflag > 0 && all(abs(newf - foundLocal) > 5e-11)
            foundLocal = [foundLocal;newf];
            if foundLocal(end) < 5e-10
                stop = true;
            end
        end
end