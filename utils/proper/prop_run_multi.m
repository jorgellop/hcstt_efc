function [img, pixx] = prop_run_multi(flnm, wl, nx, varargin)
% function [img, pixx] = prop_run_multi(flnm, wl, nx, prm)
% Execute multiple instances of a Proper prescription in parallel.
% Accepts multiple wavelengths and/or multiple parameters.
% img  = array of output images returned by the Proper prescription.
% pixx = sampling of img (m / pixel)
%        It is the responsibility of the prescription to return this
%        value (which is returned from prop_end).  If either "wl" or
%        "prm" are arrays, then "pixx" will be an array with each
%        element corresponding to the respective entry in that/those
%        arrays.
% flnm = filename (excluding extension) of Matlab Proper prescription
% wl   = scalar or vector containing wavelength(s) (um)
%        If a vector, each wavelength will be run in parallel and the
%        resulting field will be three-dimensional, with the third
%        dimension being the field at the corresponding wavelength.
%        NOTE: Unlike prop_run, this entry CANNOT be a file pointing to
%        a table of wavelengths and weights.
% nx   = size of computational grid (same for all instances)
% prm  = array of structures containing input parameters
%        (passed to Proper prescription)(optional).
% Note: if both wl and prm are specified as arrays, they must have the
% same number of elements.

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2014 Mar     jek  created idl routine
% 2015 Apr 06  gmg  Matlab translation
% The DO_TABLE, PRINT_INTENSITY, and VERBOSE options to prop_run
% are not supported by this routine.

% Default parameters
  if nargin < 1
    error('Proper:PROP_RUN_MULTI', 'Need filename of Proper prescription');
  end
  if nargin < 2
    error('Proper:PROP_RUN_MULTI', 'Need wavelength(s)');
  end
  if nargin < 3
    error('Proper:PROP_RUN_MULTI', 'Need size of computational grid');
  end
  nwl  = size(wl, 2);           % number of wavelength values
  if nargin > 3
    prm  = varargin{1};
    npc  = size(prm, 2);        % number of parameter columns
    if(( nwl > 1 ) & ( npc > 1 ) & ( npc ~= nwl ))
      error('Proper:PROP_RUN_MULTI', ...
        'Number of parameters ~= number of wavelengths');
    end
  else
    npc  = 1;                   % number of parameter columns
  end

  funh = str2func(flnm);        % function handle for Proper prescription
  if nwl > npc
    nj   = nwl;                 % number of jobs = number of wavelengths
  else
    nj   = npc;                 % number of jobs = number of parameters
  end

% Create vector of wavelength values
  if nwl > 1
    wlv  = wl;
  else
    wlv(1 : nj) = wl;
  end

% Create vector of parameter values
  if nargin > 3
    if npc > 1
      prmv = prm;
    else
      prmv(1 : nj) = prm;
    end
  end

% Find the number of logical processors available
  if ispc == 0                  % if not a Microsoft Windows OS
    [stat, cpus] = system('getconf _NPROCESSORS_ONLN');
% Note: This works in Linux and Mac OSX
%   [stat, cpus] = system('nproc');
% only works in Linux
%   [stat, cpus] = system('sysctl hw.logicalcpu');
% only works in Mac OSX
  else                          % Microsoft Windows OS
    [stat, cpus] = system('echo %NUMBER_OF_PROCESSORS%');
  end

  ncpu = str2num(cpus);         % number of available logical CPUs
  fprintf(1, '   available logical CPUs:%6d\n', ncpu);
  if nj > ncpu
    error('Proper:PROP_RUN_MULTI', ...
      'Requested cpus %6d > available cpus %6d\n', nj, ncpu);
  end

% Force Matlab 'local' cluster to recognize all available logical CPUs
  clst = parcluster('local');
  clst.NumWorkers = ncpu;
  saveProfile(clst);

% Set time scales for different operating systems
  if ismac == 1                 % Mac OSX
    tscl = 1.0e-9;              % time scale (sec)
  elseif isunix == 1            % Unix or Linux OS
    tscl = 1.0e-6;              % time scale (sec)
  else                          % Microsoft Windows OS
    tscl = 1.0e-6;              % time scale (sec)
  end

  tstc = tic;                   % start time to create pool
  pool = parpool(nj);           % create pool of processors
  tdec = toc(tstc);             % delta time to create pool
  fprintf(1, ' start time  delta time\n');
  fprintf(1, '    (sec)       (sec)\n');
  fprintf(1, '            %12.3f  create pool\n', tdec);

% Run Proper prescription (flnm) in parallel on nj logical processors
  img  = NaN(nx, nx, nj);       % one frame array
  pixx = zeros(nj, 1);          % spacing between points in x (m)
  tdei = zeros(nj, 1);          % delta time for each parallel instance
  tsti = zeros(nj, 1, 'uint64');% start time for each parallel instance
  tsts = tic;                   % start time for parallel section
  for ij = 1 : nj
    tsti(ij) = tic;
    if ( nargin < 4 )           % no parameters
      parf(ij) = parfeval(pool, funh, 2, wlv(ij), nx);
    else                        % 1 or more parameters
      parf(ij) = parfeval(pool, funh, 2, wlv(ij), nx, prmv(ij));
    end
  end
  for ij = 1 : nj
    [ijc, imgj, pixj] = fetchNext(parf);
    img(:, :, ijc) = imgj(:, :);
    pixx(ijc) = pixj;
    tdei(ijc) = toc(tsti(ijc));
    fprintf(1, '%12.6f%12.3f%6d\n', ...
      tscl * cast(tsti(ijc) - tsts, 'double'), tdei(ijc), ijc);
  end
  tdes = toc(tsts);             % delta time for parallel section
  fprintf(1, '            %12.3f  parallel section\n', tdes);

  tstd = tic;                   % start time to delete pool
  delete(pool);                 % delete pool of processors
  tded = toc(tstd);             % delta time to delete pool
  fprintf(1, '            %12.3f  delete pool\n', tded);

end                             % function prop_run_multi.m
