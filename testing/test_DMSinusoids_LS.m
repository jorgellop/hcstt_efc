% test_DMSinusoids_LS.m
%
%
%
% Jorge Llop - Mar 11, 2019
% close all;
apRad = 128/2;
Nact = 30;
pixxact = apRad*2/Nact;

actxcyc = pixxact*4;
N = 1024;                           %DM axis base value (12x12)
Y = repmat([-N/2:(N/2-1)],N,1);  %Matrix with ascending rows
X = Y.';                          %Matrix with ascending columns
RHO = sqrt(X.^2 + Y.^2);            %Matrix of element dist. from cent.
TTA = atan2(Y,X);                   %Matrix of element angle from cent.

dm_actuators_mat0 = 0.5*1*(1+cos(2*pi.*RHO.*cos(TTA-0)/actxcyc+0));

dm_actuators_mat = zeros(N);
dm_actuators_mat(RHO<apRad) = dm_actuators_mat0(RHO<apRad)*1;

aux = zeros(N);
aux(RHO<apRad) = 1;
wf = aux.*exp(1i*2*pi*dm_actuators_mat);
wf = fftshift(wf);
im = fft2(wf);
im = fftshift(im);
imagesc(abs(im).^2)
axis image
% imagesc(dm_actuators_mat)


