function [bm] = prop_begin(diam, wl, np, bdf, wfa, dx)
% [bm] = prop_begin(diam, wl, np, bdf, wfa, dx)
% Initialize variables for PROPER routines.  This routine must be called
% before any other PROPER routines in order to initialize required
% variables.
% Creates a structure with a complex array of electric field magnitudes
% = wfa (if present) or = (1.0, 0.0i) otherwise
% bm      = beam structure created by this routine
% diam    = diameter of beam (m)
% wl      = wavelength (m)
% np      = number of points in each dimension of np x np array
% bdf     = beam diameter fraction (optional) (not used if dx is specified)
% wfa     = complex np x np array of E-field values (optional)
% dx      = grid sampling (m / pixel) (optional)

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine 
% 2014 Apr 18  gmg  Matlab translation
% 2015 Apr 02  gmg  Added optional wfa input

% beam structure:
% bm.diam     = initial beam diameter (m)
% bm.dx       = grid sampling (m / pixel)
% bm.fr       = focal ratio
% bm.PropType = 'In__to_Out', 'Out_to_In_', or 'In__to_In_'
% bm.pz       = beam position z (m)
% bm.Rbeam    = beam radius of curvature (m)
% bm.RbeamInf = beam starts out with infinite curvature radius
% bm.RefSurf  = reference surface type: 'planar' or 'spheri'
% bm.TypeOld  = beam location 'In_' or 'Out' (inside or outside beam waist)
% bm.w0       = beam waist radius (m)
% bm.w0_pz    = beam waist position z (m)
% bm.wf       = initialized wavefront array
% bm.wl       = wavelength (m)
% bm.zRay     = Rayleigh distance for current beam (m)

  propcommon
  npa         = np;     % number of points across array
  OldOPD      = 0.0;
  RayFact     = 2.0;    % wavefront switches from planar to spherical
                        % at RayFact * bm.zRay

  if nargin < 1
    fprintf(1, ...
      'prop_begin: beam diam, wavelength, and number of points required\n');
    return
  elseif nargin < 2
    fprintf(1, 'prop_begin: wavelength, and number of points required\n');
    return
  elseif nargin < 3
    fprintf(1, 'prop_begin: number of points required\n');
    return
  end

  bm.diam     = diam;
  if nargin < 4
    bm.dx     = diam / 0.5 / np;        % default bdf = 0.5
  elseif nargin < 6
    bm.dx     = diam / bdf / np;
  else
    bm.dx     = dx;
  end
  bm.fr       = 1.0d9;          % default is collimated beam
  bm.PropType = 'In__to_In_';
  bm.pz       = 0.0;
  bm.Rbeam    = 0.0;
  bm.RbeamInf = 1;              % beam starts with infinite curvature radius
  bm.RefSurf  = 'PLANAR';
  bm.TypeOld  = 'In_';
  bm.w0       = diam / 2.0;
  bm.w0_pz    = 0.0;
  if nargin < 5
    bm.wf     = complex(ones(np, np), zeros(np, np));
  else
    bm.wf     = prop_shift_center(wfa, 1);
  end
  bm.wl       = wl;
  bm.zRay     = pi * bm.w0^2 / wl;

  nlist = 1500;
  if do_table
    ActionNum = 0;              % list index
    bmdl = zeros(nlist, 1);     % list of beam diameters at each lens (m)
    dzl  = zeros(nlist, 1);     % list of propagation distances (m)
    efrl =  ones(nlist, 1);     % effective focal ratios after each lens
    fll  = zeros(nlist, 1);     % list of lens focal lengths (m)
    saml = zeros(nlist, 1);     % list of sampling at each surface (m)
    snml = cellstr(char(zeros(nlist, 1)));      % list of surface names
  end

end                     % function prop_begin
