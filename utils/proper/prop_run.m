function [img, pixx] = prop_run(flnm, wl, nx, varargin)
% function [img, pixx] = prop_run(flnm, wl, nx, prm, poff, pint, tbl, qt, vrbs)
% Execute one or more instances of a Proper prescription.
% Accepts a wavelength or a file of wavelengths and weights.
% img  = array of output images returned by the Proper prescription.
% pixx = sampling of img (m / pixel)
%        It is the responsibility of the prescription to return this
%        value (which is returned from prop_end).
% flnm = filename (excluding extension) of Matlab Proper prescription
% wl   = Either the wavelength (um) or the name of a text file
%        containing a list of wavelength and weight pairs.
%        In the latter case, the prescription is run for each wavelength
%        and the results added together with the respective weights.
% nx   = size of computational grid.  Most efficient if a power of 2.
% prm  = structure containing input parameters (optional)
%        (passed to Proper prescription).
% poff = phase offset
%      = 1, then a phase offset is added as the wavefront is propagated.
%        For instance, if a wavefront is propagated over a distance of
%        1/4 wavelength, a phase offset of pi/2 radians will be added.
%        This is useful in cases where the offset between separate beams
%        that may be combined later may be important (e.g., the
%        individual arms of an interferometer).
%        By default, a phase offset is not applied.
% pint = print intensity
%      = 1, then print total intensity; default = 0
% qt   = quiet
%      = 1, then intermediate messages and surface labels will not be
%        printed; default = 0.
% tbl  = do_table
%      = 1, then prints out a table of sampling and beam size for each
%        surface; default = 0.
% vrbs = verbose
%      = 1, then informational messages will be sent; default = 0.

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine
% 2015 Jun 22  gmg  Matlab translation

  propcommon

% Default parameters
  if nargin < 1
    error('Proper:PROP_RUN', 'Need filename of Proper prescription');
  end
  if nargin < 2
    error('Proper:PROP_RUN', 'Need wavelength(s)');
  end
  if nargin < 3
    error('Proper:PROP_RUN', 'Need size of computational grid');
  end
  if nargin > 3                 % prm is given
    prm  = varargin{1};
  end
  if nargin > 4                 % poff is given
    prop_phase_offset = varargin{2};
  else
    prop_phase_offset = 0;
  end
  if nargin > 5                 % pint is given
    print_total_intensity = varargin{3};
  else
    print_total_intensity = 0;
  end
  if nargin > 6                 % qt is given
    print_it = ~varargin{4};
  else
    print_it = 1;
  end
  if nargin > 7                 % tbl is given
    do_table = varargin{5};
  else
    do_table = 0;
  end
  if nargin > 8                 % vrbs is given
    prop_verbose = varargin{6};
  else
    prop_verbose = 0;
  end

  funh = str2func(flnm);        % function handle for Proper prescription
  tic;
  if ischar(wl)                 % then wl is the name of a text file
    dlmt = ' ';                 % delimiter used in file
    nhdr = 1;                   % number of header lines in file
    wla  = importdata(wl, dlmt, nhdr);  % read wavelengths and weights
    wla.data(:, 1) = wla.data(:, 1) * 1.0e-6;
    nwl  = size(wla.data, 1);   % number of wavelengths in file
    img = zeros(nx, nx);
    for iwl = 1 : nwl
      if print_it               % print wavelength and throughput
        fprintf(1, 'Lambda(m) = %12.6e   Throughput = %8.6f\n', ...
          wla.data(iwl, 1), wla.data(iwl, 2));
      end
      if nargin > 3             % prm is given
        [imgt, pixx] = funh(wla.data(iwl), nx, prm);
      else
        [imgt, pixx] = funh(wla.data(iwl), nx);
      end
      img  = img + imgt * wla.data(iwl, 2);
    end
  else                          % wl is a single wavelength
    wl   = wl * 1.0e-6;         % convert wl to meters
    if nargin > 3               % prm is given
      [img, pixx] = funh(wl, nx, prm);
    else
      [img, pixx] = funh(wl, nx);
    end
  end
  if do_table
    prop_table;
  end
  tdet = toc;

  if print_it                   % print elapsed time
    fprintf(1, 'Total elapsed time (seconds) =      %12.6f\n', tdet);
  end
end                             % function prop_run
