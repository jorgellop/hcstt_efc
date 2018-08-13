function bm = prop_circular_aperture(bm, ra, cx, cy, norm)
% [bm] = prop_circular_aperture(bm, ra, cx, cy, norm)
% Multiplies the wavefront in bm by a circular aperture
% bm   = beam structure (input and output)
% ra   = radius of aperture (meters unless norm = 1)
% cx   = center coordinate x (meters unless norm = 1) (optional)
% cy   = center coordinate y (meters unless norm = 1) (optional)
% norm = 1 for ra, cx, cy normalized to beam radius

% 2005 Feb     jek  created idl routine
% 2014 Jun 09  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed prop_shift_center call to allow odd size arrays

  if nargin < 3
    cx   = 0.0;
    cy   = 0.0;
  end
  if nargin < 5
    norm = 0;
  end

  bm.wf = bm.wf.*prop_shift_center(prop_ellipse(bm,ra,ra,cx,cy,0.0,norm),1);
end                     % function prop_circular_aperture
