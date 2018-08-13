function w = generateTukeyWindow( Nwindow, RHO, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nlut = round(10*Nwindow);
p = linspace(-Nwindow/2,Nwindow/2,Nlut);
lut = tukeywin(Nlut,alpha);

w = interp1(p,lut,RHO,'linear',0);
end

