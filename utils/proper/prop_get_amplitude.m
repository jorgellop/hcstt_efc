function [amp] = prop_get_amplitude(bm)
% [amp] = prop_get_amplitude(bm)
% Return the amplitude of the wavefront array
% bm   = beam structure
% amp  = amplitude of the wavefront array

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation

  amp  = prop_shift_center(abs(bm.wf));
end                     % function prop_get_amplitude
