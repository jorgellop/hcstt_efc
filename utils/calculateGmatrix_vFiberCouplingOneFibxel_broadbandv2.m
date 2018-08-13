function [ G ] = calculateGmatrix_vFiberCouplingOneFibxel_broadbandv2( bmb4dm,usePix, surf_DM1_init, Q, Nact, ac_spac, poke_amp, infl, lambda0, N , info)
%calculateGmatrix calculates the G matrix for EFC
%THIS IS FOR THE USE OF COMBINED G
%   Detailed explanation goes here
    
    normIc = info.normIc;
%     fibermode0 = info.fibermode0;
%     Q_fib = info.Q_fib;
    
    lam_arr = info.lam_arr;
    if(usePix)
        fiberDiam_pix = info.fiberDiam_pix4pix;
        fibermode0 = info.fibermode04pix;
        factor_fiber_thruput = info.factor_fiber_thruput4pix;
        use_fiber = false;
        normal_EFC = true;
    else
        fiberDiam_pix =  info.fiberDiam_pix;
        fibermode0 = info.fibermode0;
        factor_fiber_thruput = info.factor_fiber_thruput;
        use_fiber = true;
        normal_EFC = false;
    end
    
    numOfWavelengths = info.numOfWavelengths;
    if(normal_EFC)
%         G = [];
%         for II = 1:numOfWavelengths
            G = [zeros(nnz(Q),Nact^2)];
%         end
    else
        G = zeros(numOfWavelengths,Nact^2); % One fibxel
    end
	count = 0; 
    
    
    % propagate to image plane using compact model and initial DM shape
%     wf_aux = fftshift(bmb4dm).*fftshift(exp(1i*4*pi*surf_DM1_init/(800e-9)));
%     bm0 = fftshift(fft2(ifftshift(wf_aux)));
    bm0 = prescription_DM1toImage_compact_vFiberCoupling_broadband_combG( bmb4dm,usePix, surf_DM1_init, true, info);

    % Poke each actuator 
    for ix = 1:Nact
        for iy = 1:Nact
            count = count + 1;
            disp(['Poking actuator ',num2str(count),'/',num2str(Nact^2)]);
            DM_strokes = zeros(N);

            xpos = round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
            ypos = round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
    
            DM_strokes(xpos,ypos) = poke_amp;
            delta_surf_DM = conv2(DM_strokes,infl,'same');            
            
            % Propagate to image plane using compact model using poked DM
            % surface height
%             wf_aux = fftshift(bmb4dm).*fftshift(exp(1i*4*pi*(surf_DM1_init+delta_surf_DM)/(800e-9)));
%             bm = fftshift(fft2(ifftshift(wf_aux)));
            bm = prescription_DM1toImage_compact_vFiberCoupling_broadband_combG( bmb4dm,usePix, surf_DM1_init+delta_surf_DM, true, info);
%             bm(Q_fib) = bm(Q_fib).*fibermode0(:);
            % Extract the field in the dark hole and vectorize (FP_col = Eab in Amir's SPIE paper)
            FP0_col = [];
            FP_col = [];
            FP_colSMF = [];
            FP_colSMF0 = [];
            for index = 1:numel(lam_arr)
                lam = lam_arr(index);
                fibermode0JJ = fibermode0(:,:,index);
                fibermode0JJ_col = fibermode0JJ(Q(:,:,index) == 1);

                bm0_lam = bm0(:,:,index);
                bm_lam = bm(:,:,index);
                
                FP0 = bm0_lam;
                FP = bm_lam;

                FP0_col = [FP0_col;(lambda0/lam)*FP0(Q(:,:,index)==1)/sqrt(normIc)];
                FP_col = [FP_col;(lambda0/lam)*FP(Q(:,:,index)==1)/sqrt(normIc)];
%                 FP_colSMF = [FP_colSMF;(lambda0/lam)*FP(Q(:,:,index)==1).*fibermode0JJ_col/sqrt(normIc)];
%                 FP_colSMF0 = [FP_colSMF0;(lambda0/lam)*FP0(Q(:,:,index)==1).*fibermode0JJ_col/sqrt(normIc)];
                if(use_fiber)
                    G(index,count) = sum(sum((lambda0/lam)*FP(Q(:,:,index)==1).*fibermode0JJ_col/sqrt(normIc)))-sum(sum((lambda0/lam)*FP0(Q(:,:,index)==1).*fibermode0JJ_col/sqrt(normIc))); % The column of the G matrix is the difference in the field 
                end
            end
            if(normal_EFC)
                G(:,count) = FP_col-FP0_col; % The column of the G matrix is the difference in the field 
%                 G(:,count) = FP_colSMF-FP_colSMF0; % The column of the G matrix is the difference in the field 
            end
%             else
%                 for JJ = 1:numOfWavelengths
%                     if(use_fiber)
%                         G(JJ,count) = sum(sum(FP_colSMF(1+fiberDiam_pix(JJ,1)^2*(JJ-1):fiberDiam_pix(JJ,1)^2*(JJ))))-sum(sum(FP_colSMF0(1+fiberDiam_pix(JJ,1)^2*(JJ-1):fiberDiam_pix(JJ,1)^2*(JJ)))); % The column of the G matrix is the difference in the field 
%                     else
%     %                     G(1,count) = sum(sum(FP_col*factor_fiber_thruput))-sum(sum(FP0_col*factor_fiber_thruput)); % The column of the G matrix is the difference in the field
%                         disp('Does not work')
%                         return;
%                     end
%                 end
%             end
        end
    end
    if(info.useGPU)
        G = gather(G);
    end
end

