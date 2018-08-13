function bm = prop_rectangular_obscuration(bm, sx, sy, cx, cy, rot, norm)
% bm = prop_rectangular_obscuration(bm, sx, sy, cx, cy, rot, norm)
% Multiplies the wavefront in bm by a rectangular obscuration
% bm   = beam structure (input and output)
% sx   = size along x of obscuration (meters unless norm = 1)
% sy   = size along y of obscuration (meters unless norm = 1)
% cx   = center coordinate x (meters unless norm = 1) (optional)
% cy   = center coordinate y (meters unless norm = 1) (optional)
% rot  = counter-clockwise rotation of rectangle about its center (deg)(opt)
% norm = 1 for sx, sy, cx, cy normalized to beam diameter

% 2005 Feb     jek  created idl routine
% 2016 Jan 14  gmg  Matlab translation

  if nargin < 4
    cx   = 0.0;
    cy   = 0.0;
  end
  if nargin < 6
    rot  = 0;
  end
  if nargin < 7
    norm = 0;
  end

  bm.wf = bm.wf .* ...
          prop_shift_center(prop_rectangle(bm,sx,sy,cx,cy,rot,1.0,norm),1);
end                     % function prop_rectangular_obscuration
