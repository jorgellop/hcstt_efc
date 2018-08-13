function [ Q ] = generateDarkHoleMask( darkHoleShape, IWA, OWA, numOfWavelengths, lambdas, lambda0, lambdaOverD, N )
%generateDarkHoleMask Summary of this function goes here
%   Detailed explanation goes here

    % Initialize variables 
    [X,Y] = meshgrid(-N/2:N/2-1); 
    [THETA,RHO] = cart2pol(X,Y);

    % Defines the dark hole
    Q = zeros(N,N,numOfWavelengths);
    for index = 1:numOfWavelengths
        lam = lambdas(index); 
        scal = lambda0/lam;
        if(strcmp(darkHoleShape,'full'))
            Q(:,:,index) = and(RHO >= IWA*lambdaOverD*scal, RHO <= OWA*lambdaOverD*scal); % Full annulus
        elseif(strcmp(darkHoleShape,'half'))
            Q(:,:,index) = and(RHO >= IWA*lambdaOverD*scal, RHO <= OWA*lambdaOverD*scal);
            Q(:,:,index) = and(Q(:,:,index), X > 0);% Half annulus
        elseif(strcmp(darkHoleShape,'60deg_half'))
            Q(:,:,index) = and(RHO >= IWA*lambdaOverD*scal, RHO <= OWA*lambdaOverD*scal);
            Q(:,:,index) = and(Q(:,:,index), X > 0);
            Q(:,:,index) = and(Q(:,:,index), abs(THETA) < pi/6); % 60deg keystone about x-axis
        elseif(strcmp(darkHoleShape,'30deg_half'))
            Q(:,:,index) = and(RHO >= IWA*lambdaOverD*scal, RHO <= OWA*lambdaOverD*scal);
            Q(:,:,index) = and(Q(:,:,index), X > 0);
            Q(:,:,index) = and(Q(:,:,index), abs(THETA) < pi/12); % 30deg keystone about x-axis
        elseif(strcmp(darkHoleShape,'Dshape'))
            Q(:,:,index) = and(X >= IWA*lambdaOverD*scal, RHO <= OWA*lambdaOverD*scal);% D shaped dark hole
        elseif(strcmp(darkHoleShape,'SmallCirc'))

            cent = (IWA+OWA)/2*scal*lambdaOverD;
            Q(:,:,index) = sqrt((X-cent).^2+Y.^2) <= 1*lambdaOverD*scal;
        elseif(strcmp(darkHoleShape,'ET'))

            im = round(rgb2gray(imread('ETmovieposter.jpg'))/255);
            croprows = 274:368;
            cropcols = 54:225;
            im = im(croprows,cropcols);
            im = bwareaopen(im, 500);
            im = imresize(im,scal);
            
            [rows0,cols0] = size(im);
        
            Q(:,:,index) = padarray_centered(flipud(im),rows0,cols0,N);
        else
            error('Undefined dark hole');
        end
    end
end

