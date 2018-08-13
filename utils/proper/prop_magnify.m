function aro = prop_magnify(ari, magn, csrv, sox, quik)
% function aro = prop_magnify(ari, magn, csrv, sox, quik)
% Resample input array using damped sinc or cubic interpolation
% Calls either an external C routine to do a damped sinc interpolation
% or uses Matlab's interp2 routine to do a cubic interpolation.
% Outputs:
% aro  = 2D output array
% Required inputs:
% ari  = input array to be magnified
% magn = magnification (e.g., 0.5 = shrink array by factor of 2)
% Optional inputs:
% csrv = 0 = Image amplitudes are not scaled. (default)
%      = 1 = The intensity in the array will be conserved.
%          = The array is assumed to be an electric field, so the
%            interpolated result will be divided by the magnification.
%      = 2 = The intensity of the array will be conserved.
%          = The array is assumed to be an intensity, so the
%            interpolated result will be divided by the magnification^2.
% sox  = size out of output array x (pixels)
%        default: the input array dimensions are multiplied by the
%        magnification
% quik = 1 = Use the external damped sinc interpolation routine (default)
%      = 3 = Use the Matlab interp2 cubic interpolation routine

% 2005 Feb     jek  created idl routine
% 2014 Aug 19  gmg  Matlab translation
% 2015 Mar 18  gmg  Added call to interp2

  if nargin < 2
    error('Proper:PROP_MAGNIFY', ...
          'Must specify input array and magnification.\n');
  end

  if nargin < 4
    sox  = size(ari, 2) * magn; % size of output array x
    soy  = size(ari, 1) * magn; % size of output array y
  else
    soy  = sox;
  end

  if nargin < 5
    quik = 1;                   % use external C damped sinc interpolation
  end

  if(quik == 3)                 % use Matlab interp2 cubic interpolation
                                % works for non-square arrays
    [ny, nx] = size(ari);       % size of input array
    cx   = fix(nx / 2) + 1;     % pixel coordinates of ari center x
    cy   = fix(ny / 2) + 1;     % pixel coordinates of ari center y
    sx   = fix(sox / 2) + 1;    % coordinates of aro center x
    sy   = fix(soy / 2) + 1;    % coordinates of aro center y
    [aox, aoy] = meshgrid(([1:sox]-sx) / magn + cx, ([1:soy]-sy) / magn + cy);
    if isreal(ari) == 0         % ari is complex
      aror = zeros(soy, sox);
      aror = interp2(real(ari), aox, aoy, 'linear');
      aroi = zeros(soy, sox);
      aroi = interp2(imag(ari), aox, aoy, 'linear');
      aro  = complex(aror, aroi);
    else                        % ari is not complex
      aro  = zeros(soy, sox);
      aro  = interp2(ari, aox, aoy, 'linear');
    end
  else                          % use external C damped sinc interpolation
                                % assumes square arrays
    if isreal(ari) == 0         % ari is complex
      aror = zeros(sox, sox);
      aror = prop_szoom(real(ari), size(ari, 2), sox, magn);
      aroi = zeros(sox, sox);
      aroi = prop_szoom(imag(ari), size(ari, 2), sox, magn);
      aro  = complex(aror, aroi);
    else                        % ari is not complex
      aro  = zeros(sox, sox);
      aro  = prop_szoom(ari, size(ari, 2), sox, magn);
    end
  end                           % if(quik == 3)

% Conserve intensity?
  if nargin < 3 | csrv == 0     % do not scale array amplitudes
  elseif csrv == 1              % scale array amplitudes by magn
    aro  = aro / magn;
  elseif csrv == 2              % scale amplitudes by magn^2
    aro  = aro / magn^2;
  end

end                     % function prop_magnify
