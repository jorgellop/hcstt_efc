function [ fun ] = fct2min_sensingEFwFiber(x,matfit,fibermode0,Eab_fib_re)

aux = x+matfit;
aux = aux(:)';
fun = aux*fibermode0(:)-Eab_fib_re;

end