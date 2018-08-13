function bm = prescription_DM1toImage_compact_vFiberCoupling_broadband( bm_BB,surf_DM1, useFPM, info )
%prescription_DM1toImage_compact Calculates image of a star
%   Detailed explanation goes here

%!!!!!!!! ONLY WORKS FOR ONE DM (IN PUPIL) !!!!!!!!!!!!!!

    apRad = info.apRad; 
    lambdaOverD = info.lambdaOverD;
    RHO = info.RHO;
    N = info.N;
    lam_arr = info.lam_arr;
    
    if(info.useApodizer)
        APOD = info.APOD;
    end
    FPM = info.FPM;
    LPM = info.LPM;

    useGPU = info.useGPU;
    for index = 1:numel(lam_arr)
        bm_lam = fftshift(bm_BB(:,:,index));
        
        lam = lam_arr(index);

        if(useGPU)
            bm_lam = gpuArray(bm_lam);
        end
        
        % apply DM1
        bm_lam = bm_lam.*fftshift(exp(1i*4*pi*surf_DM1/lam));
        
        if(info.useApodizer)
            bm_lam = bm_lam.*fftshift(APOD);
        end
        
%         figure(1);
%         imagesc(angle(fftshift(bm_lam.wf)));
%         axis image; 
%         colorbar; 
%         drawnow;
        
        if(useFPM)
            LP = vortexCoronagraph_Pup2Pup( fftshift(bm_lam), FPM, 0, apRad, lambdaOverD, RHO, N, 'dft', 'forward', useGPU );
        else
            LP = vortexCoronagraph_Pup2Pup( fftshift(bm_lam), ones(N), 0, apRad, lambdaOverD, RHO, N, 'fft', 'forward', useGPU );
        end

%         figure(2);
%         imagesc(abs(LP.*LPM));
%         axis image; 
%         colorbar; 
%         drawnow;
        
        bm_lam = circshift(rot90(myfft2(LP.*LPM),2),[1,1]);
        bm_BB(:,:,index) = bm_lam;
%         figure(3);
%         imagesc(abs(bm_lam.wf));
%         axis image; 
%         colorbar;    
%         drawnow;
    end
    bm = bm_BB;
end

