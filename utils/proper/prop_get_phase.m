function [ph] = prop_get_phase(bm)
% [ph] = prop_get_phase(bm)
% Return the phase of the wavefront array
% bm   = beam structure
% ph   = phase of the wavefront array

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation

  ph   = prop_shift_center(angle(bm.wf));
end                     % function prop_get_phase
