function [ PUPIL, scale_factor ] = loadPupilFile( fullpath, scale_factor, apRad, N, findRadius)
%loadPupilFile Summary of this function goes here
%   Detailed explanation goes here

    PUPIL = fitsread(fullpath);
    [rows0,cols0] = size(PUPIL);
    if(findRadius)
        [rowslist,colslist] = find(round(PUPIL));
        apRadin = max(sqrt((rowslist-rows0/2-1).^2 + (colslist-cols0/2-1).^2))
        clear rowslist colslist
        scale_factor = apRad/apRadin
        
    end
    
    PUPIL = imresize(PUPIL,scale_factor,'lanczos3');
    [rows0,cols0] = size(PUPIL);
    PUPIL = padarray_centered(PUPIL,rows0,cols0,N);
    
end

