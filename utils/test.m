clear; close all; 

addpath('segMirrorFunctions');

N = 1024; % Size of NxN computational grid 
apRad = 500; % aperture radius: (flat to flat diameter)/2 in samples 
apDia = 2*apRad; % flat to flat diameter in samples 
gapWidth = 0.02/12*apDia;% gap width in samples 
numRings = 3;% number of rings of hexagonal segments (for Keck, numRings=3)

%% Generate pupil function

hexMirror.apDia = apDia; % flat to flat aperture diameter (samples)
hexMirror.gapWidth = gapWidth; % samples
hexMirror.numRings = numRings;% Number of rings in hexagonally segmented mirror 
hexMirror.N = N;

[ PUPIL ] = hexSegMirror_getSupport( hexMirror );

figure;
imagesc(PUPIL);
colormap(gray);
colorbar; 
axis image; axis xy; 
title('pupil support');

%% Generate pupil field

numSegments = hexSegMirror_numSegments( numRings );

hexMirror.pistons = rand(1,numSegments);% Vector of pistons per segment (waves)
hexMirror.tiltxs = rand(1,numSegments); % Vector of x-tilts per segment (waves/apDia)
hexMirror.tiltys = rand(1,numSegments);% Vector of y-tilts per segment (waves/apDia)

[ PUPIL ] = hexSegMirror_getField( hexMirror );

figure;
ax1=subplot(1,2,1);
imagesc(abs(PUPIL));
colormap(ax1,gray(256));
colorbar; 
axis image; axis xy; 
title('amplitude');

ax2=subplot(1,2,2);
imagesc(angle(PUPIL));
colormap(ax2,hsv(256));caxis([-pi pi]);
colorbar; 
axis image; axis xy; 
title('phase');
