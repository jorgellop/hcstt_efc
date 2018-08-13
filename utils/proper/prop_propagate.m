function bm = prop_propagate(bm, dz, snm, to_plane)
% [bm] = prop_propagate(bm, dz, snm, to_plane)
% Determine which propagator to use to propagate the current wavefront
% by a specified distance and do it.
% bm       = beam structure (input and output)
% dz       = distance to propagate wavefront (m)
% snm      = name of surface to which to propagate (string) (optional)
% to_plane = if true, propagation is to a plane

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine
% 2014 May 09  gmg  Matlab translation

  propcommon

  if print_it
    if nargin > 2
      fprintf(1, 'Propagating to %s\n', snm);
    else
      fprintf(1, 'Propagating\n');
    end
  end

  [bm, dzw]   = prop_select_propagator(bm, dz);
  z2   = bm.pz + dz;

  if nargin > 3 & to_plane
    bm.PropType = [bm.PropType(1 : 7) 'In_'];
  end

  if prop_verbose
    fprintf(1, '  Propagate: prop_type:   %10s\n', bm.PropType);
  end

  if     strcmp(bm.PropType, 'In__to_In_')
    bm   = prop_ptp(bm, dz                 );

  elseif strcmp(bm.PropType, 'In__to_Out')
    bm   = prop_ptp(bm, bm.w0_pz - bm.pz   );
    bm   = prop_wts(bm, z2       - bm.w0_pz);

  elseif strcmp(bm.PropType, 'Out_to_In_')
    bm   = prop_stw(bm, bm.w0_pz - bm.pz   );
    bm   = prop_ptp(bm, z2       - bm.w0_pz);

  elseif strcmp(bm.PropType, 'Out_to_Out')
    bm   = prop_stw(bm, bm.w0_pz - bm.pz   );
    bm   = prop_wts(bm, z2       - bm.w0_pz);
  end

  if print_total_intensity
    tint = sum(sum(abs(bm.wf).^2));
    if nargin > 2
      fprintf(1, 'Total intensity at surface %s: %10.3e\n', snm, tint);
    else
      fprintf(1, 'Total intensity:          %10.3e\n', tint);
    end
  end

% Save stuff for layout plots

  if do_table
    ActionNum = ActionNum + 1;  % list index
    bmdl(ActionNum) = 2 * prop_get_beamradius(bm);
%   bmdl = list of beam diameters at each lens (m)
    dzl(ActionNum)  = dz;       % list of propagation distances (m)
    saml(ActionNum) = bm.dx;    % list of sampling at each surface (m)
    if nargin > 2
      snml{ActionNum} = snm;    % list of surface names
    else
      snml{ActionNum} = '(Surface)';
    end
  end
end                     % function prop_propagate
