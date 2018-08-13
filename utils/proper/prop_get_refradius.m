function [rr] = prop_get_refradius(bm)
% [rr] = prop_get_refradius(bm)
% Return the radius of the reference sphere to which the current
% wavefront's phase is reference.  If the reference surface is planar
% (near field), then the radius will be 0.0.  Assuming that forward
% propagation occurs from left to right, a negative radius indicates
% that the center of the reference surface is to the right of the
% current position (e.g., the beam is converging).
% The reference radius is defined to be the distance from the pilot beam
% waist to the current position.
% bm   = beam structure
% rr   = radius of curvature of sphere to which the current wavefront's
%        phase is referenced (m)

% 2006 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation

  rr   = bm.pz - bm.w0_pz;
end                     % function prop_get_refradius
