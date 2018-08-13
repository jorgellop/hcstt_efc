% function [bm, wf] = ProperUse(wl, np, diam, flob, flep, fley)
% Runs the simple telescope example from Proper in Matlab
% bm   = beam structure created by this routine
% wf   = wavefront array at end
% wl   = wavelength (m)
% np   = number of points in each dimension of np x np array
% diam = diameter of beam (m)
% flob = focal length of objective lens of telescope (m)
% flep = focal length of eyepiece of telescope (m)
% fley = focal length of eye (m)

% 2014 Jun 10  gmg

  diam = 0.060;                 % diameter of the objective (m)
  flob = 15.0 * diam;           % focal length of the objective (m)
  flep = 0.021;                 % focal length of the eyepiece (m)
  fley = 0.022;                 % focal length of the eye (m)
  wl   = 500e-9;                % wavelength (m)
  np   = 1024;                  % number of points

  epd  = flep / (1.0 - flep / (flob + flep));   % exit pupil distance (m)

  bm = prop_begin(diam, wl, np);

  bm = prop_circular_aperture(bm, diam / 2.0);

  bm = prop_define_entrance(bm);

  bm = prop_lens(bm, flob);

  bm = prop_propagate(bm, flob + flep);

  bm = prop_lens(bm, flep);

  bm = prop_propagate(bm, epd);

  bm = prop_lens(bm, fley);

  bm = prop_propagate(bm, fley);

  figure(1); imagesc(abs(bm.wf)); axis xy; axis equal; colorbar; colormap(gray)

  figure(2); imagesc(angle(bm.wf)); axis xy; axis equal; colorbar

  wf = prop_end(bm);

  figure(3); imagesc(log10(wf)); axis xy; axis equal; colorbar; colormap(gray)

  figure(4); plot(wf(np / 2 + 1, 1 : np));
