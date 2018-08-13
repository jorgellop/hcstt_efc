function [bm] = prop_add_phase(bm, perr)
% [bm] = prop_add_phase(bm, perr)
% Add a phase error map or value to the current wavefront array.
% The phase error array is assumed to be at the same sampling as the
% current wavefront.  Note that this is wavefront, not surface, error.
% bm   = wavefront structure (input and output)
% perr = scalar or 2D image containing phase error in meters

% 2005 Feb     jek  created idl routine
% 2014 May 08  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed prop_shift_center call to allow odd size arrays

  if size(perr, 1) == 1
    bm.wf = bm.wf .* exp(2.0 * pi * i * perr / bm.wl);
  else
    bm.wf = bm.wf .* exp(2.0 * pi * i * prop_shift_center(perr, 1) / bm.wl);
  end
end                     % function prop_add_phase
