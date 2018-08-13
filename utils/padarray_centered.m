function [ output_image ] = padarray_centered( image, rows, cols, N)
%padarray_centered Pads image such that it is centered in NxN grid
%   Detailed explanation goes here

% pad arrays to size NxN
if(and(mod(rows,2)==0,mod(cols,2)==0))
    image = padarray(image,[(N-rows)/2 (N-cols)/2]);
elseif(mod(rows,2)==0)
    image = padarray(image,[(N-rows)/2 0]);
    image = padarray(image,[0 floor((N-cols)/2)],0,'pre');
    image = padarray(image,[0 floor((N-cols)/2)+1],0,'post');
elseif(mod(cols,2)==0)
    image = padarray(image,[0 (N-cols)/2]);
    image = padarray(image,[floor((N-rows)/2) 0],0,'pre');
    image = padarray(image,[floor((N-rows)/2)+1 0],0,'post');
else
    image = padarray(image,[0 floor((N-cols)/2)],0,'pre');
    image = padarray(image,[0 floor((N-cols)/2)+1],0,'post');    
    image = padarray(image,[floor((N-rows)/2) 0],0,'pre');
    image = padarray(image,[floor((N-rows)/2)+1 0],0,'post');
end
output_image = image; 
end

