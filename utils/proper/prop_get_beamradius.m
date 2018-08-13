function br = prop_get_beamradius(bm)
% [br] = prop_get_beamradius(bm)
% This routine returns the half-width-at-half-max of the Gaussian beam
% that is used to keep track of the beam size.  The Gaussian beam starts
% off with a beam diameter equal to the entrance pupil diameter.
% This ignores widening of the beam due to aberrations.
% bm   = beam structure
% br   = beam radius (m)

% 2005 Feb     jek  created idl routine
% 2014 May 09  gmg  Matlab translation

  br   = bm.w0 * sqrt(1.0 + (bm.wl * (bm.pz - bm.w0_pz) / (pi * bm.w0^2))^2);
end                     % function prop_get_beamradius
