function [ c, ceq ] = nonlcon_sensingEFwFiber(x,im_fib_pixres,scale,matfit,fibermode0)

c = [];

sz = size(x);
n = sz(1)/scale;

count = 0;
for II = 1:n
    for JJ = 1:n
        count = count + 1;
        eab_fit = 0;
        fibparam = 0;
        for KK = 1:scale
            for LL = 1:scale
%                 eab_fit = eab_fit + (x(scale*(II-1)+KK,scale*(JJ-1)+LL) ...
%                      + matfit(scale*(II-1)+KK,scale*(JJ-1)+LL)) ...
%                      * fibermode0(scale*(II-1)+KK,scale*(JJ-1)+LL);
%                 fibparam = fibparam + fibermode0(scale*(II-1)+KK,scale*(JJ-1)+LL);
                eab_fit = eab_fit + x(scale*(II-1)+KK,scale*(JJ-1)+LL) ...
                     + matfit(scale*(II-1)+KK,scale*(JJ-1)+LL);
            end
        end
%         ceq(II,JJ) = eab_fit - im_fib_pixres(II,JJ)*fibparam/scale^2;
        ceq(II,JJ) = eab_fit - im_fib_pixres(II,JJ);
    end
end