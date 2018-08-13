function [wfa, samp, sampsec] = prop_end(bm, idx, noabs)
% [wfa, samp, sampsec] = prop_end(bm, idx, noabs)
% Set variables needed to properly conclude a propagation run
% wfa     = wavefront array (output)
% samp    = wavefront array sampling distance (m) (optional output)
% sampsec = wavefront array sampling distance (arcsec) (optional output)
% bm      = beam structure (input)
% idx     = if nonzero, returns the idx by idx central part of wavefront
% noabs   = returns intensity if 0 (default), complex wavefront if 1

% 2005 Feb     jek  created idl routine
% 2014 May 07  gmg  Matlab translation

  samp    = prop_get_sampling(bm);
  sampsec = prop_get_sampling_arcsec(bm);

  if nargin == 3 & noabs == 1           % keep complex E-field values
    wfa   = prop_shift_center(bm.wf);
  else                                  % calculate intensities
    wfa   = prop_shift_center(abs(bm.wf).^2);
  end

  if nargin >= 2 & idx ~= 0             % select central portion of array
    idxn =  ceil((idx - 1) / 2);
    idxp = floor((idx - 1) / 2);
    i2   = floor(size(bm.wf) / 2) + 1;  % indices of center of array
    wfa  = wfa(i2(1) - idxn : i2(1) + idxp, i2(2) - idxn : i2(2) + idxp);
  end
end                     % function prop_end
