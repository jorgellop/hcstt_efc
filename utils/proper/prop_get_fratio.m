function [fr] = prop_get_fratio(bm)
% [fr] = prop_get_fratio(bm)
% Return the current wavefront focal ratio 
% This routine computes the current beam's focal ratio by dividing
% the current distance to focus by the current beam diameter.
% bm   = beam structure
% fr   = focal ratio

% 2005 Feb     jek  created idl routine
% 2014 May 07  gmg  Matlab translation

  fr   = bm.fr;
end                     % function prop_get_fratio
