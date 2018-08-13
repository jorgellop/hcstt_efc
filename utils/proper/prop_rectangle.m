function apm = prop_rectangle(bm, sx, sy, cx, cy, rot, dark, norm)
% function apm = prop_rectangle(bm, sx, sy, cx, cy, rot, dark, norm)
% Return an image containing an antialiased, filled rectangle
% apm  = aperture mask containing antialiased filled rectangle
% bm   = beam structure
% sx   = size along x (meters unless norm = 1, then fraction of beam diameter)
% sy   = size along y (meters unless norm = 1, then fraction of beam diameter)
% cx   = center of rect. relative to wf center x (m unless norm = 1)(optional)
% cy   = center of rect. relative to wf center y (m unless norm = 1)(optional)
% rot  = counter-clockwise rotation of rectangle about its center (deg)(opt)
% dark = 1 = draw a dark rect. (0 inside, 1 outside)(default is opposite way)
% norm = 1 = sizes and center coordinates are normalized to beam diameter

% 2005 Feb     jek  created idl routine
% 2015 Nov 11  gmg  Matlab translation

  dx   = prop_get_sampling(bm);         % spacing between points in x (m)
  dy   = prop_get_sampling(bm);         % spacing between points in y (m)
  magn = 11;                    % subsampling factor at edges; must be odd
  nx   = size(bm.wf, 2);                % number of pixels in x
  ny   = size(bm.wf, 1);                % number of pixels in y
  prx  = prop_get_beamradius(bm) / dx;  % beam radius x in pixels
  pry  = prop_get_beamradius(bm) / dy;  % beam radius y in pixels

  if nargin < 4
    cpx  = floor(nx / 2 + 1);           % center in x (pixels)
    cpy  = floor(ny / 2 + 1);           % center in y (pixels)
    radx = 0.5 * sx / dx;               % radius in x (pixels)
    rady = 0.5 * sy / dy;               % radius in y (pixels)
  elseif nargin > 7 & norm == 1
    cpx  = floor(nx / 2 + 1) + cx * prx;
    cpy  = floor(ny / 2 + 1) + cy * pry;
    radx = 0.5 * sx * prx;              % radius in x (pixels)
    rady = 0.5 * sy * pry;              % radius in y (pixels)
  else
    cpx  = floor(nx / 2 + 1) + cx / dx;
    cpy  = floor(ny / 2 + 1) + cy / dy;
    radx = 0.5 * sx / dx;               % radius in x (pixels)
    rady = 0.5 * sy / dy;               % radius in y (pixels)
  end

  if nargin < 6
    rot = 0.0;                          % rotation (degrees)
  end

  if nargin < 7
    dark = 0;
  end

% Calculate rectangle vertex x and y coordinates before rotation and shift
  vx0  = [-radx, -radx,  radx,  radx];  % (pixels)
  vy0  = [-rady,  rady,  rady, -rady];  % (pixels)
  nvrt = 4;                             % number of vertices

% Calculate rectangle vertex x and y coordinates after rotation and shift
  vx   = vx0 * cosd(rot) - vy0 * sind(rot) + cpx;       % (pixels)
  vy   = vx0 * sind(rot) + vy0 * cosd(rot) + cpy;       % (pixels)

  apm  = zeros(ny, nx);                 % initialize mask to zeros

  left = find(vy == min(vy));
  left = left(find(vx(left) == min(vx(left))));
  left = left(1);                       % index of lower left corner (1 - 4)
  if left ~= nvrt
    ltnt = left + 1;                    % index of next left corner (1 - 4)
  else
    ltnt = 1;                           % index of next left corner (1 - 4)
  end
  rght = left;
  if rght ~= 1
    rtnt = rght - 1;                    % index of next right corner (1 - 4)
  else
    rtnt = nvrt;                        % index of next right corner (1 - 4)
  end

  ivyn = round(min(vy));                % index of vertex minimum y
  if ivyn < 2
    ivyn = 2;
  end
  ivyx = round(max(vy));                % index of vertex maximum y
  if ivyx > (ny - 1)
    ivyx = ny - 1;
  end

  for ipy = ivyn : ivyx                 % index of pixel y
    for suby = 0 : (magn - 1)
% Calculate y coordinate of current edge pixel (pixels)
      py   = ipy - 0.5 + (0.5 + suby) / magn;

      if py < vy(left)
        continue;
      end
      if py > max(vy)
        break;
      end
      if py >= vy(ltnt)
        left = ltnt;
        if left ~= nvrt
          ltnt = left + 1;
        else
          ltnt = 1;
        end
      end

      if py >= vy(rtnt)
        rght = rtnt;
        if rght ~= 1
          rtnt = rght - 1;
        else
          rtnt = nvrt;
        end
      end

      dylt = vy(ltnt) - vy(left);
      if dylt ~= 0
        dxlt = vx(ltnt) - vx(left);
        xlt  = dxlt / dylt * (py - vy(left)) + vx(left);
      else
        xlt  = vx(left);
      end

      dyrt = vy(rtnt) - vy(rght);
      if dyrt ~= 0
        dxrt = vx(rtnt) - vx(rght);
        xrt  = dxrt / dyrt * (py - vy(rght)) + vx(rght);
      else
        xrt  = vx(rght);
      end

      xltp = round(xlt);
      xrtp = round(xrt);

      if xltp ~= xrtp
        if xltp >= 1
          apm(ipy, xltp) = apm(ipy, xltp) + magn * ((xltp + 0.5) - xlt);
        end
        if xrtp <= nx
          apm(ipy, xrtp) = apm(ipy, xrtp) + magn * (xrt - (xrtp - 0.5));
        end
        if (xrtp - xltp) > 1
          xlp1 = max(xltp + 1,  1);
          xrp1 = min(xrtp - 1, nx);
          apm(ipy, xlp1 : xrp1) = apm(ipy, xlp1 : xrp1) + magn;
        end
      else
        if ((xltp >= 1) & (xltp <= nx))
          apm(ipy, xltp) = apm(ipy, xltp) + magn * (xrt - xlt);
        end
      end                               % if xltp ~= xrtp
    end                                 % for suby = 0 : (magn - 1)
  end                                   % for ipy = ivyn : ivyx

  apm  = apm / magn^2;
  if dark == 1
    apm  = 1.0 - apm;
  end
