N = 2^10; 

% Initialize variables 
[X,Y] = meshgrid(-N/2:N/2-1); 

numOfWavelengths = 3; 
scals = [0.9,1.0,1.1];

darkHoleShape = 'ET';

    % Defines the dark hole
    Q = zeros(N,N,numOfWavelengths);
    for index = 1:numOfWavelengths

            scal = scals(index);
            im = round(rgb2gray(imread('ETmovieposter.jpg'))/255);
            croprows = 274:368;
            cropcols = 54:225;
            im = im(croprows,cropcols);
            im = bwareaopen(im, 500);
            im = imresize(im,scal);
            [rows0,cols0] = size(im);
          
            Q(:,:,index) = padarray_centered(flipud(im),rows0,cols0,N);
    end


    figure;
    imagesc(Q(:,:,1));
    colorbar;
    colormap(gray);
    axis image
    
    figure;
    imagesc(Q(:,:,2));
    colorbar;
    colormap(gray);
    axis image
    
    figure;
    imagesc(Q(:,:,3));
    colorbar;
    colormap(gray);
    axis image

% im = rgb2gray(imread('ETmovieposter.jpg'));
% 
% croprows = 274:368;
% cropcols = 54:225;
% 
% im = round(im(croprows,cropcols)/255);
% 
% [rows0,cols0] = size(im);
% 
% Q = padarray_centered(im,rows0,cols0,N);
% 
% imagesc(Q)
% 
% % Q = zeros(N);
% % Q = insertText(Q,position,text_str,'FontSize',18,'BoxColor',...
% %     'k','BoxOpacity',0,'TextColor','white','AnchorPoint','center);
% % 
% % imshow(Q)
% % 
% % % I = imread('peppers.png');
% % % J = insertText(I, [100 315 ], 'Peppers are good for you!');
% % % imshow(J);