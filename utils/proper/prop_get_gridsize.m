function [nx, ny] = prop_get_gridsize(bm)
% [nx, ny] = prop_get_gridsize(bm)
% Return the dimensions of the current wavefront
% bm   = beam structure
% nx   = grid size in x
% ny   = grid size in y

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation

  nx   = size(bm.wf, 2);
  ny   = size(bm.wf, 1);
end                     % function prop_get_gridsize
