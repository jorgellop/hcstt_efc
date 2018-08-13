function bm = prop_stw(bm, dz)
% [bm] = prop_stw(bm, dz)
% Propagate a wavefront from a spherical reference surface that is
% outside the Rayleigh limit from its focus to a planar reference
% surface that is inside.  Used by prop_propagate.
% Intended only for use by prop_propagate.  Not a user-callable routine.
% bm   = beam structure (input and output)
% dz   = distance to propagate wavefront (m)

% 2005 Feb     jek  created idl routine
% 2014 May 12  gmg  Matlab translation

  if abs(dz) <= 1e-12
    return
  end

  global prop_verbose
  if prop_verbose
    fprintf(1, '  STW: dz:                %10.3f\n', dz);
  end

  if ~strcmp(bm.RefSurf, 'SPHERI')
    fprintf(1, '  STW: Input reference surface not spherical; using ptp.\n');
    bm = prop_ptp(bm, dz);
    return
  end
  if nargin < 2
    dz   = bm.w0_pz - bm.pz;
  end
  bm.pz = bm.pz + dz;
  bm.dx = bm.wl * abs(dz) / size(bm.wf, 1) / bm.dx;

  nx   = size(bm.wf, 2);
  ny   = size(bm.wf, 1);
  srn  = sqrt(nx * ny);

  if dz > 0.0           % use forward transform
    bm.wf =  fft2(bm.wf);
% Note that  fft2 includes a normalization of 1
    bm.wf = bm.wf / srn;
  else                  % use inverse transform
    bm.wf = ifft2(bm.wf);
% Note that ifft2 includes a normalization of N^-2
    bm.wf = bm.wf * srn;
  end

  bm   = prop_qphase(bm, dz);

  global prop_phase_offset
  if prop_phase_offset
    bm.wf = bm.wf * exp(i * 2 * pi * dz / bm.wl);
  end

  if prop_verbose
    fprintf(1, '  STW: z:                 %10.3f'  , bm.pz);
    fprintf(1, '       dx:                %10.3e\n', bm.dx);
  end

  bm.RefSurf = 'PLANAR';
end                     % function prop_stw
