%{
DM Writing Function: Write Sinusoid
- Places a Sin Function on Mirror Surface
- Sin calculated using DE_DMMapSin.m function
- Utilizes DE_DMArrayToVect to shape map for write
*** ASSUMES MIRROR CONNECTION ALREADY PRESENT; DOES NOT CLOSE CONNECTION
*** Defaults to DE_DMMapSin output setting 7: only returns heigh
        If desired, hnm can be plotted with imagesc. Refer to DE_DMMapSin
 
******************************************************
- Arguments:
    h0          = Max poke height in nm
    q           = angle of sinusoid
    x0          = actuators per cycle
    alp         = phase delay
    drv_info    = DM info from OPEN_mutliDM
- Returns:
    hnm         = Surface map as matrix in nm
    hV          = Vector of voltage percentages for writing
******************************************************

Compiled By:    Daniel Echeverri
Last Modified:  08/04/2016
%}

function im = hcstt_TakeModelImage(u,FPM,info)

apRad = info.apRad;
normalize = info.normalize;
N = 1024;
Nact = info.Nact;

if info.posDM_x ~= 0
    ac_spac = info.ac_spac;
else
    ac_spac = round(2*apRad/Nact);
end
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac );

posDM_x = info.posDM_x;
posDM_y = info.posDM_y;

count = 0; 
DM1_strokes = zeros(N,N);
for ix = 1:Nact
    for iy = 1:Nact
        count = count + 1;
        if(posDM_x(1) == 0 && posDM_y(1) == 0)
            xpos = round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
            ypos = round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
        else
            xpos = round(posDM_x(ix));%
            ypos = round(posDM_y(iy));%%
        end
        
        DM1_strokes(xpos,ypos) = u(count) + DM1_strokes(xpos,ypos);
    end
end
surf_DM1 = conv2(DM1_strokes,infl,'same');

[X,Y] = meshgrid(-N/2:N/2-1); 
[THETA,RHO] = cart2pol(X,Y);
lambdaOverD = N/apRad/2; % lambda/D (samples) 

info.apRad = apRad; 
info.lambdaOverD = lambdaOverD;
info.RHO = RHO;
info.N = N;
info.lam_arr = 650e-9;
info.useGPU = false;
info.FPM = exp(1i*8*THETA);
info.LPM = exp(-(RHO/(0.85*apRad)).^1000);
info.useApodizer = false;

% wf1 = complex(ones(N, N), zeros(N, N));
wf1 = exp( -(RHO/(apRad)).^100 );
wf2 = prescription_DM1toImage_compact_vFiberCoupling_broadband( wf1, surf_DM1, FPM, info);

im = abs(wf2).^2;
if(normalize)
    im = im/max(im(:));
end
end
