function bm = prop_rectangular_aperture(bm, sx, sy, cx, cy, rot, norm, rnd)
% bm = prop_rectangular_aperture(bm, sx, sy, cx, cy, rot, norm, rnd)
% Multiplies the wavefront in bm by a rectangular aperture
% bm   = beam structure (input and output)
% sx   = size along x of aperture (meters unless norm = 1)
% sy   = size along y of aperture (meters unless norm = 1)
% cx   = center coordinate x (meters unless norm = 1) (optional) (default = 0)
% cy   = center coordinate y (meters unless norm = 1) (optional) (default = 0)
% rot  = counter-clockwise rotation of rectangle about its center (deg)(opt)
% norm = 1 for sx, sy, cx, cy normalized to beam diameter
% rnd  = radius of rounded corners (optional) (m)
% Note: cannot specify both rotation and rounded corners

% 2005 Feb     jek  created idl routine
% 2016 Jan 19  gmg  Matlab translation

  if nargin == 8 & rot ~= 0.0
    error('Proper:PROP_RECTANGULAR_APERTURE','Cannot specify both ROT and RND');
  end
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

  if nargin < 8
    bm.wf = bm.wf .* ...
            prop_shift_center(prop_rectangle(bm,sx,sy,cx,cy,rot,0.0,norm),1);
  end
end                     % function prop_rectangular_aperture
