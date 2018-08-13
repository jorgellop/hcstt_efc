function bm = prop_elliptical_aperture(bm, rx, ry, cx, cy, norm)
%        bm = prop_elliptical_aperture(bm, rx, ry, cx, cy, norm)
% Multiplies the wavefront in bm by an elliptical clear aperture
% bm   = beam structure (input and output)
% rx   = radius of aperture x (meters unless norm = 1)
% ry   = radius of aperture y (meters unless norm = 1)
% cx   = center coordinate x (meters unless norm = 1) (optional)
% cy   = center coordinate y (meters unless norm = 1) (optional)
% norm = 1 for rx, ry, cx, cy normalized to beam radius (optional)

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine
% 2016 Apr 14  gmg  Matlab translation

  if nargin < 4
    cx   = 0.0;
    cy   = 0.0;
  end
  if nargin < 6
    norm = 0;
  end

  bm.wf = bm.wf.*prop_shift_center(prop_ellipse(bm,rx,ry,cx,cy,0.0,norm),1);
end                     % function prop_elliptical_aperture
