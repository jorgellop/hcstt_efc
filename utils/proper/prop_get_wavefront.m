function [wf] = prop_get_wavefront(bm)
% [wf] = prop_get_wavefront(bm)
% Return the complex-valued wavefront array
% bm   = beam structure
% wf   = complex-valued wavefront array

% 2005 Feb     jek  created idl routine
% 2014 Jul 29  gmg  Matlab translation

  wf   = prop_shift_center(bm.wf);
end                     % function prop_get_wavefront
