function bm = prop_qphase(bm, cu)
% [bm] = prop_qphase(bm, cu)
% This routine applies a quadratic, radially-dependent phase factor
% caused by a wavefront curvature to the current wavefront.
% It is used by prop_stw and prop_wts.
% Intended only for internal use by prop_* routines.
% bm   = beam structure (input and output)
% cu   = phase curvature radius (m)

% 2005 Feb     jek  created idl routine
% 2013 Oct     jek  speed up phase calculation
% 2014 May 12  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed fftshift to ifftshift to allow for odd size arrays

  if cu == 0.0
    return
  end

  bm.wf = bm.wf .* exp(i * pi * ifftshift(prop_radius2(bm)) / cu / bm.wl);
end                     % function prop_qphase
