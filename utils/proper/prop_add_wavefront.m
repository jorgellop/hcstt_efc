function [bm] = prop_add_wavefront(bm, wfad)
% [bm] = prop_add_wavefront(bm, wfad)
% Add a wavefront to the current wavefront array.
% The wavefront array is assumed to be at the same sampling as the
% current wavefront.  Note that this is wavefront, not surface, error.
% bm   = wavefront structure (input and output)
% wfad = scalar or 2D image containing the value or wavefront to add.
%        NOTE: All responsibility is on the user to ensure that the
%        two fields have the same sampling and reference phase curvature.
%        NO CHECKING is done by this routine.

% 2011 Mar     jek  created idl routine
% 2016 Apr 05  gmg  Matlab translation

  if size(wfad, 1) == 1
    bm.wf = bm.wf + wfad;
  else
    bm.wf = bm.wf + prop_shift_center(wfad, 1);
  end
end                     % function prop_add_wavefront
