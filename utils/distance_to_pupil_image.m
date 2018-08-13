function [ z2 ] = distance_to_pupil_image( z1, f1, f2 )
%distance_to_pupil_image Summary of this function goes here
%   Detailed explanation goes here

z2 = f2*(f1*(f1+f2)-f2*z1)/f1^2;

end

