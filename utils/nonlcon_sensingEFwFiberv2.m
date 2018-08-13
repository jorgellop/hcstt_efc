function [ c, ceq ] = nonlcon_sensingEFwFiberv2(x,matfit,fibermode0,Eab_fib_re)

c = [];
aux = x+matfit;
aux = aux(:)';
ceq = aux*fibermode0(:)-Eab_fib_re;

end