function [mapo] = prop_resamplemap(bm, mapi, pxsc, cpx, cpy, shx, shy)
% [mapo] = prop_resamplemap(bm, mapi, pxsc, cpx, cpy, shx, shy))
% Interpolate input map using cubic convolution onto grid with same size
% and sampling as the current wavefront array.  Optionally shift the map.
% mapo = output map
% bm   = beam structure
% mapi = aberration map to be resampled
% pxsc = spacing of "mapi" (m)
% cpx  = pixel coordinate of mapi center x (0,0 is center of 1st pixel)
% cpy  = pixel coordinate of mapi center y (0,0 is center of 1st pixel)
% shx  = amount to shift map x (m) (optional)
% shy  = amount to shift map y (m) (optional)

% 2005 Feb     jek  created idl routine
% 2014 Jul 29  gmg  Matlab translation

  if nargin < 7
    shx  = 0.0;
    shy  = 0.0;
  end

  [ny, nx] = size(bm.wf);       % number of points in wavefront array
  o1x  = ([1 : nx] - fix(nx / 2) - 1)  * bm.dx / pxsc;
  o1x  = o1x + cpx - shx / pxsc;
  o1y  = ([1 : ny] - fix(ny / 2) - 1)' * bm.dx / pxsc;
  o1y  = o1y + cpy - shy / pxsc;

  mapo = interp2(mapi, o1x, o1y, 'cubic', 0.0);

end                     % function prop_resamplemap
