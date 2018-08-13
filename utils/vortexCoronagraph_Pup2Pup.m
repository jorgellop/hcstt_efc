    function OUT = vortexCoronagraph_Pup2Pup( IN, FPM, R0, apRad, lambdaOverD, RHO, N, algo, operation, useGPU )
%vortexCoronagraph Summary of this function goes here
%   Detailed explanation goes here

    if(useGPU)
        IN = gpuArray(IN);
    end
    if( strcmp(algo,'fft') && ~strcmp(operation,'adj') )
        
        EP = IN;
        OUT = myifft2(myfft2(EP).*FPM);
        
	elseif( strcmp(algo,'fft') && strcmp(operation,'adj') )
        
        LP = IN;
        OUT = myifft2(myfft2(LP).*conj(FPM));
        
    elseif( strcmp(algo,'dft') )
        
        if(R0 == 0)
            cut_rad1 = 10.1735/pi;
            cut_rad2 = 32.1897/pi;
        else
            obj = @(x) besselj(1,x) - R0*besselj(1,R0*x);
            cut_rad1 = fzero(obj,10)/pi;
            cut_rad2 = fzero(obj,32)/pi;
        end

        windowKnee = 1-cut_rad1/cut_rad2;

        D = 2*apRad;
        
        % !!!!!!!!!!! I have forced the maximum number of samples to work with proper
        
        NA = N;
        crop = N/2-NA/2+1:N/2+NA/2;
        
        
        % DFT vectors 
        x = ((0:NA-1)-NA/2)/D;
        u1 = ((0:N-1)-N/2)/lambdaOverD;
        u2 = ((0:N-1)-N/2)*2*cut_rad2/N;
        
        windowMASK1 = generateTukeyWindow( 2*cut_rad2*lambdaOverD, RHO, windowKnee ) ;
        windowMASK2 = generateTukeyWindow( N, RHO, windowKnee ) ;
        
        if(useGPU)
            x = gpuArray(x);
            u1 = gpuArray(u1);
            u2 = gpuArray(u2);
            windowMASK1 = gpuArray(windowMASK1);
            windowMASK2 = gpuArray(windowMASK2);
        end
        if(~strcmp(operation,'adj'))

            EP = IN;
            EP = EP(crop,crop);

            %%%%%%% Large scale DFT

            FP1 = (N/lambdaOverD)/(D*N)*exp(-1i*2*pi*u1'*x)*EP*exp(-1i*2*pi*x'*u1); 
            LP1 = (N/lambdaOverD)/(D*N)*exp(1i*2*pi*x'*u1)*(FP1.*FPM.*(1-windowMASK1))*exp(1i*2*pi*u1'*x);

            %%%%%%% Fine sampled DFT

            FP2 = 2*cut_rad2/(D*N)*exp(-1i*2*pi*u2'*x)*EP*exp(-1i*2*pi*x'*u2); 
            LP2 = 2*cut_rad2/(D*N)*exp(1i*2*pi*x'*u2)*(FP2.*FPM.*windowMASK2)*exp(1i*2*pi*u2'*x);

            OUT = padarray_centered(LP1+LP2,NA,NA,N);
            %disp('Propagating through vortex with forward DFT.');
            
        elseif(strcmp(operation,'adj') )
        
            LP = IN(crop,crop);

            %%%%%%% Large scale DFT

            FP1 = (N/lambdaOverD)/(D*N)*exp(-1i*2*pi*u1'*x)*LP*exp(-1i*2*pi*x'*u1); 
            EP1 = (N/lambdaOverD)/(D*N)*exp(1i*2*pi*x'*u1)*(FP1.*conj(FPM).*(1-windowMASK1))*exp(1i*2*pi*u1'*x);


            %%%%%%% Fine sampled DFT

            FP2 = 2*cut_rad2/(D*N)*exp(-1i*2*pi*u2'*x)*LP*exp(-1i*2*pi*x'*u2); 
            EP2 = 2*cut_rad2/(D*N)*exp(1i*2*pi*x'*u2)*(FP2.*conj(FPM).*windowMASK2)*exp(1i*2*pi*u2'*x);

            OUT = padarray_centered(EP1+EP2,NA,NA,N);
        end
    else
        error('Error. \nChoose algo = fft or dft. \nChoose forward or adj operation.')
    end
    if(useGPU)
        OUT = gather(OUT);
    end
end

