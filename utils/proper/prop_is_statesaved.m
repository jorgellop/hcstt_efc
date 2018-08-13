function stat = prop_is_statesaved(bm)
%        stat = prop_is_statesaved(bm)
% stat = status = 1 if a state exists, 0 otherwise
% bm   = current beam structure
% Determine if a previously saved state exists for the current wavelength.

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine
% 2016 Apr 06  gmg  Matlab translation

  propcommon

  stat = 0;
  if numel(save_state_lam) ~= 0
    for il = 1 : size(save_state_lam, 1)
% Does a state exist for the current wavelength?
      if save_state_lam(il) == bm.wl
        stat = 1;
      end
    end
  end
end                     % function prop_is_statesaved
