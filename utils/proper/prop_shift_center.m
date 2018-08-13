function arrs = prop_shift_center(arr, inv)
% [arrs] = prop_shift_center(arr, inv)
% Shift an array from the origin (1, 1) to the center (default)
% or from the center to the origin (inv = 1)
% arr  = complex or real array (input)
% arrs = complex or real array shifted (output)

% 2005 Feb     jek  created idl routine
% 2014 May 08  gmg  Matlab translation
% 2014 Sep 02  gmg  Added inverse shift to correctly handle odd size arrays

  if nargin == 2 & inv == 1
    arrs = ifftshift(arr);
  else
    arrs =  fftshift(arr);
  end
end                     % function prop_shift_center
