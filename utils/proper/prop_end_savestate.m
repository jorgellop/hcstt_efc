function prop_end_savestate
%        prop_end_savestate
% Terminate the current save state system.
% This deletes the files created by prop_state/prop_savestate.
% This routine should be called only once during a session.

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine
% 2016 Apr 06  gmg  Matlab translation

  propcommon

  save_state = 0;
  save_state_lam = 0;

  system(['rm *' statefile]);
end                     % function prop_end_savestate
