function [ wf_fib ] = wfOverFibxel( wf, Q_fib)

numlam = numel(wf(1,1,:));
wf_fib = zeros(nnz(Q_fib(:,:,1)),numlam);
for II = 1:numlam
    wf_II = wf(:,:,II);
    wf_fib(:,II) = wf_II(Q_fib(:,:,II)==1);
end
end