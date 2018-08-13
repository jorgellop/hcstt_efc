function [dmcz, dmsf] = prop_fit_dm(dmz, infk)
% [dmcz, dmsf] = prop_fit_dm(dmz, infk)
% Determine deformable mirror actuator piston values
% that generate a desired DM surface, accounting for
% the effect of the actuator influence function.
% dmcz = DM actuator positions that create the desired surface
%        when the influence function is included (2D image)
% dmsf = dmcz array convolved with infk
% dmz  = DM surface to match (2D array, with each element representing
%        the desired surface amplitude for the corresponding actuator)
% infk = influence function kernel (2D image sampled at actuator spacing)
% Intended for use by the prop_dm routine and not for general users

% 2005 Feb     jek  created idl routine
% 2014 Jun 25  gmg  Matlab translation

  e0   = 1.0e6;                         % initial assumed rms error

  last = dmz;                           % last good dm
  dmsf = conv2(dmcz, infk, 'same');     % 2D convolution
  diff = dmz - dmsf;
  erms = sqrt(sum(sum(diff^2)));        % root mean square error

  while erms < e0 & (e0 - erms) / erms > 0.01
    last = dmcz;                        % last good dm
    e0   = erms;                        % previous iteration rms error
    dmcz = dmcz + diff;
    dmsf = conv2(dmcz, infk, 'same');
    diff = dmz - dmsf;
    if max(max(abs(diff))) < 1.0e-15
      break
    end
    erms = sqrt(sum(sum(diff^2)));
  end

% If fit diverged, use the result from the previous iteration
  if erms > e0
    dmcz = last;
  end

end                             % function prop_fit_dm
