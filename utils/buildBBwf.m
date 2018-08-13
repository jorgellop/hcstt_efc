% Builds the wavefront with multiple arrays accounting for the different
% wavelengths
function [ bm_BB ] = buildBBwf(planet,lambda0,lam_arr,N,RHO,apRad,error_map,x_fib,lambdaOverD)
wf1_noerrors = complex(ones(N, N), zeros(N, N));
for II = 1:numel(lam_arr)
    scal = lambda0/lam_arr(II);
    wf1_noerrorsII = wf1_noerrors;
%     wf1_noerrorsII(RHO > apRad * scal) = 0;
    wf1_noerrorsII =  exp( -(RHO/(apRad*scal)).^100 );
%     figure(2);
%     imagesc(abs(wf1_noerrorsII(N/2-80:N/2+80,N/2-80:N/2+80)).^2)
%     axis image
%     colorbar
%     drawnow

    if(planet)
        [X,Y] = meshgrid(-N/2:N/2-1); 
        wf1_noerrorsII = wf1_noerrorsII.*exp(-1i * 2 * pi *x_fib*lambdaOverD*X/N );
    end

    lam = lam_arr(II);
    wf1 = fftshift(wf1_noerrorsII).*fftshift(exp(1i*4*pi*error_map/lam));
    wf1 = fftshift(wf1);
%     figure(2);
%     imagesc(abs(wf1_noerrorsII(N/2-80:N/2+80,N/2-80:N/2+80)).^2)
%     axis image
%     colorbar
%     drawnow

    bm_BB(:,:,II) = wf1;
end
end