function [bm, dzw] = prop_select_propagator(bm, dz)
% [bm, dzw] = prop_select_propagator(bm, dz)
% Used by prop_propagate to decide in which propagation regime
% the next surface will be (to decide which propagation method to use
% (spherical-to-planar, planar-to-spherical, or planar-to-planar)).
% bm   = beam structure (input and output)
% dz   = distance to new position z from current position z (m)
% dzw  = distance to new waist position from new position z (m)

% 2005 Feb     jek  created idl routine
% 2014 May 09  gmg  Matlab translation

  global RayFact

  dzw  = bm.w0_pz - bm.pz;
  newz = bm.pz + dz;

  if abs(bm.w0_pz - newz) < RayFact * bm.zRay
    TypeNew  = 'In_';
  else
    TypeNew  = 'Out';
  end

  bm.PropType = [bm.TypeOld '_to_' TypeNew];
  bm.TypeOld  = TypeNew;
end                     % function prop_select_propagator
