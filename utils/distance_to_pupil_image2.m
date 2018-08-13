function [ z1 ] = distance_to_pupil_image2( z2, f1, f2 )
%distance_to_pupil_image Summary of this function goes here
%   Detailed explanation goes here

z1 = (f1*(f1*f2 + f2^2 - f1*z2))/f2^2;

end

