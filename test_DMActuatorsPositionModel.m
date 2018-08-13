clear all;
% close all;

apRad = 2^6;

N = 1028;
Nact = 12;
ac_spac = round(2*apRad/Nact);
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac );

count = 0; 
DM1_strokes = zeros(N,N);

u = zeros(Nact,Nact);
act = [10,10];
u(act(1),act(2)) = 0.002;

for ix = 1:Nact
    for iy = 1:Nact
        count = count + 1;
        xpos = round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
        ypos = round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
        if [ix,iy] == act
            [xpos,ypos]
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
info.FPM = exp(1i*6*THETA);
info.LPM = ones(N,N);
info.useApodizer = false;

% wf1 = complex(ones(N, N), zeros(N, N));
wf1 = exp( -(RHO/(apRad)).^100 );
% wf2 = prescription_DM1toImage_compact_vFiberCoupling_broadband( wf1, surf_DM1, true, info);
        bm_lam = fftshift(wf1);
        lam = 650e-9;
        % apply DM1
        bm_lam = bm_lam.*fftshift(exp(1i*4*pi*surf_DM1/lam));
%         LP = vortexCoronagraph_Pup2Pup( fftshift(bm_lam), ones(N), 0, apRad, lambdaOverD, RHO, N, 'fft', 'forward', false );

        wf2 = myfft2(fftshift(bm_lam));
        
        wf2(:,1:N/2+3) = 0;
        
        wf3_pup = myifft2(wf2);
        
        wf4 = circshift(rot90(myfft2(wf3_pup),2),[1,1]);
        wf5_pup = myifft2(wf4);
im1 = abs(wf5_pup).^2;
% im1 = imag(wf3_pup);
im = im1/max(im1(:));

imagesc(im);
ll = N/2-100;
ul = N/2+100;
axis image
axis([ll ul ll ul])
colorbar
