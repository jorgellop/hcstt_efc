function [ G ] = calculateGmatrix( bmb4dm, surf_DM1_init, Q, Nact, ac_spac, poke_amp, infl, lambda0, N , info)
%calculateGmatrix calculates the G matrix for EFC
%   Detailed explanation goes here
    
    info.show_plots = false;
    normIc = info.normIc; % Normalization for compact model 

    G = zeros(nnz(Q),Nact^2);
	count = 0; 

    % propagate to image plane using compact model and initial DM shape
    bm0 = prescription_DM1toImage_compact( bmb4dm, surf_DM1_init, true, info);

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

%                 figure(101);
%                 imagesc(info.xvals/info.apRad,info.yvals/info.apRad,delta_surf_DM.*(info.RHO<info.apRad));
%                 axis image; 
%                 colorbar; caxis([0 poke_amp]);
%                 axis([-1 1 -1 1]);
%                 drawnow;
            
            
            % Propagate to image plane using compact model using poked DM
            % surface height
            bm = prescription_DM1toImage_compact( bmb4dm, surf_DM1_init+delta_surf_DM, true, info);

            % Extract the field in the dark hole and vectorize (FP_col = Eab in Amir's SPIE paper)
            FP0_col = [];
            FP_col = [];
            for index = 1:numel(bm)
                bm0_lam = bm0(index);
                bm_lam = bm(index);
                
                FP0 = bm0_lam.wf;
                FP = bm_lam.wf;
                
%                 figure(102);
%                 imagesc(info.xvals,info.yvals,abs(FP-FP0).^2/normIc);
%                 axis image; 
%                 colorbar; 
%                 axis([-100 100 -100 100]);
%                 caxis([0 1e-5])
%                 drawnow;

                FP0_col = [FP0_col;(lambda0/bm0_lam.wl)*FP0(Q(:,:,index)==1)/sqrt(normIc)];
                FP_col = [FP_col;(lambda0/bm_lam.wl)*FP(Q(:,:,index)==1)/sqrt(normIc)];
            end
            G(:,count) = FP_col-FP0_col; % The column of the G matrix is the difference in the field 
        end
    end
    if(info.useGPU)
        G = gather(G);
    end
end

