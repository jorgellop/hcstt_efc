function bm = prop_define_entrance(bm)
% [bm] = prop_define_entrance(bm)
% Normalize the wavefront amplitude for a total power of 1.
% bm   = beam structure (input and output)

% 2005 Feb     jek  created idl routine
% 2014 Jun 10  gmg  Matlab translation

  propcommon

  pwr  = sum(sum(abs(bm.wf).^2));
  bm.wf = bm.wf / sqrt(pwr);
end                     % function prop_define entrance
