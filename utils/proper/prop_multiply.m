function [bm] = prop_multiply(bm, mult)
% function [bm] = prop_multiply(bm, mult)
% Multiplies the wavefront array by a value or an array.
% bm   = beam structure
% mult = either a 2D array (the same size as the wavefront array) or a scalar
%        The mult array is assumed to be centered at pixel(ny/2+1, nx/2+1).

% 2005 Feb     jek  created idl routine
% 2014 Aug 29  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed prop_shift_center call to allow for odd size arrays

  [my, mx] = size(mult);
  if my * mx > 1                % mult is a 2D array
    [ny, nx] = size(bm.wf);
    if mx ~= nx | my ~= ny      % mult size does not equal wavefront size
      error('Proper:PROP_MULTIPLY', ...
      'mult array size (%d, %d) not equal to wf array size (%d, %d).\n', ...
      my, mx, ny, nx);
    end
    bm.wf = bm.wf .* prop_shift_center(mult, 1);
  else
    bm.wf = bm.wf * mult;
  end
end                     % function prop_multiply
