function bm = prop_lens(bm, fl, snm)
% [bm] = prop_lens(bm, fl, snm)
% Alter the current wavefront as a perfect lens would.
% This routine computes the phase change caused by a perfect lens
% that has a focal length specified by the user.
% A positive focal length corresponds to a convex lens or concave mirror;
% a negative focal length corresponds to a concave lens or convex mirror.
% This routine updates the new beam waist position.
% bm   = wavefront structure (input and output)
% fl   = focal length of lens (m)
% snm  = surface name (string) (optional)

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine
% 2014 May 08  gmg  Matlab translation

  propcommon
  if print_it
    if nargin < 3
      fprintf(1, 'Applying lens\n');
    else
      fprintf(1, 'Applying lens at %s\n', snm);
    end
  end

  bm.zRay = pi * bm.w0^2 / bm.wl;
  wsrf = bm.w0 * sqrt(1.0 + ((bm.pz - bm.w0_pz) / bm.zRay)^2);
% wsrf = beam radius at lens surface (m)

  dzw0 = bm.pz - bm.w0_pz;
% fprintf(1, '  Lens: dzw0:       %16.9e'  , dzw0);
% fprintf(1, '        wsrf:       %16.9e\n', wsrf);

  if dzw0 ~= 0.0        % lens is not at focus or entrance pupil
    gRbo = dzw0 + bm.zRay^2 / dzw0;
    gRboinf  = 0.0;
    if gRbo ~= fl
      gRbm = 1.0 / (1.0 / gRbo - 1.0 / fl);
      gRbminf  = 0.0;   % output beam is spherical
      if prop_verbose
        fprintf(1, '  Lens: Gaussian Rbeam old: %8.3f'  , gRbo);
        fprintf(1, '        R beam:             %8.3f\n', gRbm);
      end
    else
      gRbminf  = 1.0;   % output beam is planar
      if prop_verbose
        fprintf(1, '  Lens: Gaussian Rbeam old: %8.3f'  , gRbo);
        fprintf(1, '        R beam: infinite         \n');
      end
    end
  else                  % at focus or entrance pupil, input beam is planar
    gRboinf  = 1.0;
    gRbm = -fl;
    gRbminf  = 0.0;     % output beam is spherical
    if prop_verbose
      fprintf(1, '  Lens: Gaussian Rbeam old: infinite');
      fprintf(1, '        R beam:             %8.3f\n', gRbm);
    end
  end

  if strcmp(bm.TypeOld, 'In_') | strcmp(bm.RefSurf, 'PLANAR')
    Rbo  = 0.0;
  else
    Rbo  = dzw0;
  end

  if not(gRbminf)
    bm.w0_pz = -gRbm / (1.0 + (bm.wl * gRbm / (pi * wsrf^2))^2) + bm.pz;
    bm.w0 = wsrf / sqrt(1.0 + (pi * wsrf^2 / (bm.wl * gRbm))^2);
  else                  % output beam is planar
    bm.w0_pz = bm.pz;
    bm.w0 = wsrf;
  end
% fprintf(1, '  Lens: w0:         %16.9e'  , bm.w0);
% fprintf(1, '        w0_pz:      %16.9e\n', bm.w0_pz);

% Determine new Rayleigh distance from focus;
% if currently inside this, then output beam is planar

  bm.zRay = pi * bm.w0^2 / bm.wl;
  dzw0 = bm.pz - bm.w0_pz;
% fprintf(1, '  Lens: dzw0:       %16.9e'  , dzw0);
% fprintf(1, '        zRay:       %16.9e\n', bm.zRay);

  if abs(dzw0) < RayFact * bm.zRay
    TypeNew  = 'In_';
    Rbm  = 0.0;
  else
    TypeNew  = 'Out';
    Rbm  = dzw0;
  end

  bm.PropType = [bm.TypeOld '_to_' TypeNew];

% Apply phase changes as needed, but don't apply if the phase
% is going to be similarly altered during propagation.

  if prop_verbose
    fprintf(1, '  Lens: propagator type:  %10s\n', bm.PropType);
    if strcmp(bm.TypeOld, 'In_')
      sRbo = 'Infinite';
    else
      sRbo = sprintf('%8.3f', Rbo);
    end
    if strcmp(TypeNew, 'In_')
      sRbm = 'Infinite';
    else
      sRbm = sprintf('%8.3f', Rbm);
    end
    fprintf(1, '  Lens: Rbmold: %8s', sRbo);
    fprintf(1, '        R_beam: %8s', sRbm);
    fprintf(1, '  focal length: %8.3f\n', fl);
    fprintf(1, '  Lens: beam diam at lens:%10.3e\n', 2.0 * wsrf);
  end

  if     strcmp(bm.PropType, 'In__to_In_')
    phas = 1.0 / fl;
 %  fprintf(1, 'In__to_In_  phas:   %16.9e\n', phas);

  elseif strcmp(bm.PropType, 'In__to_Out')
    phas = 1.0 / fl + 1.0 / Rbm;
 %  fprintf(1, 'In__to_Out  phas:   %16.9e  Rbm:    %16.9e\n', phas, Rbm);

  elseif strcmp(bm.PropType, 'Out_to_In_')
    phas = 1.0 / fl - 1.0 / Rbo;
 %  fprintf(1, 'Out_to_In_  phas:   %16.9e  Rbo:    %16.9e\n', phas, Rbo);

  elseif strcmp(bm.PropType, 'Out_to_Out')
    if     Rbo == 0.0
      phas = 1.0 / fl + 1.0 / Rbm;
    elseif Rbm == 0.0
      phas = 1.0 / fl - 1.0 / Rbo;
    else
      phas = 1.0 / fl - 1.0 / Rbo + 1.0 / Rbm;
      if prop_verbose
        fprintf(1, '  Lens: 1/lens_fl:        %10.3f\n', 1.0 / fl );
        fprintf(1, '  Lens: 1/R_beam_old:     %10.3f\n', 1.0 / Rbo);
        fprintf(1, '  Lens: 1/R_beam:         %10.3f\n', 1.0 / Rbm);
        fprintf(1, '  Lens: lens_phase:       %10.3f\n', phas);
      end
    end
  end

  bm   = prop_add_phase(bm, -prop_radius2(bm) * phas / 2.0);

  if strcmp(TypeNew, 'In_')
    bm.RefSurf  = 'PLANAR';
  else
    bm.RefSurf  = 'SPHERI';
  end
  bm.TypeOld  = TypeNew;
  bm.fr       = abs(dzw0) / wsrf / 2.0;

% Save stuff for layout plots

  if do_table
    ActionNum = ActionNum + 1;  % list index
    bmdl(ActionNum) = 2 * wsrf; % list of beam diameters at each lens (m)
    efrl(ActionNum) = bm.fr;    % effective focal ratios after each lens
    fll(ActionNum)  = fl;       % list of lens focal lengths (m)
    saml(ActionNum) = bm.dx;    % list of sampling at each surface (m)
    if nargin > 2
      snml{ActionNum} = snm;    % list of surface names
    else
      snml{ActionNum} = '(Lens)';
    end
  end

  if prop_verbose
    fprintf(1, '  Lens: Rayleigh distance:  %8.2e\n', bm.zRay);
  end
end                     % function prop_lens
