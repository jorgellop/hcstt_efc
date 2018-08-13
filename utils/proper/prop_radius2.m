function r2 = prop_radius2(bm)
% [r2] = prop_radius2(bm)
% Returns a 2D array in which the value of each element corresponds to
% the square of the distance of that element from the center of the
% current wavefront.
% The center of the wavefront is set to be at the center of the array.
% bm   = beam structure
% r2   = radius^2 (m^2)

% 2005 Feb     jek  created idl routine
% 2014 May 09  gmg  Matlab translation

  nx   = size(bm.wf, 2);
  ny   = size(bm.wf, 1);
  ix2  = [-floor(nx / 2 ) : floor((nx - 1) / 2)].^2;
  iy2  = [-floor(ny / 2 ) : floor((ny - 1) / 2)].^2;
  [ax2, ay2]   = meshgrid(ix2, iy2);
  r2   = (ax2 + ay2) * bm.dx^2;
end                     % function prop_radius2
