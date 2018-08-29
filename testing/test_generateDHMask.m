lambdaOverD = 1024/142/2;
N = 1024;
x_fib = 2.5;
y_fib = 0;

throughput = 0.8119;
throughput_arr = [];
numtry  = 75;
diam_arr = linspace(0.5,3,numtry);
for JJ = 1:numtry
     Q_pixJJ = exp( -(RHO/(diam_arr(JJ)/2*lambdaOverD*scal)).^100 );
     throughput_arr = [throughput_arr; sum(sum(abs(fftshift(bm_ideal(index).wf).*Q_pixJJ).^2))/info.totalPower];
end
[THETA_fib,RHO_fib] = cart2pol(X - x_fib * lambdaOverD *scal ,Y);
diamII = interp1(throughput_arr,diam_arr,throughput);
Q(:,:,index) = exp(-(RHO_fib/(diamII/2*lambdaOverD*scal)).^100);
