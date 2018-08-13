function prop_savestate(bm)
%        prop_savestate(bm)
% bm   = beam structure
% Write out the current wavefront state to a file for the current wavelength.

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine
% 2016 Apr 06  gmg  Matlab translation

  propcommon
  if save_state ~= 0

    [ny, nx] = size(bm.wf);
    iz   = 0;
    sarr = zeros(2 * ny * nx, 1);
    for iy = 1 : ny
      for ix = 1 : nx
        iz   = iz + 1;
        sarr(iz) = real(bm.wf(iy, ix));
        iz   = iz + 1;
        sarr(iz) = imag(bm.wf(iy, ix));
      end
    end

    fno  = [num2str(bm.wl * 1.0e9) statefile];
    fid  = fopen(fno, 'w');             % file ID
    if fid == -1
      error('Proper:PROP_SAVESTATE', ...
        'Cannot open file %s.\n', fno);
    end
    fwrite(fid, pwr, 'double');         % wavefront power
    fwrite(fid, sarr, 'double');        % complex wavefront amplitude
    fwrite(fid, bm.wl, 'double');       % wavelength (m)
    fwrite(fid, bm.dx, 'double');       % grid sampling (m / pixel)
    if bm.TypeOld == 'Out'
      fwrite(fid, 'OUTSIDE', 'char');
    else
      fwrite(fid, 'INSIDE_', 'char');
    end
    if bm.RefSurf == 'spheri'
      fwrite(fid, 'SPHERI', 'char');
    else
      fwrite(fid, 'PLANAR', 'char');
    end
    fwrite(fid, bm.Rbeam, 'double');    % beam radius of curvature (m)
    fwrite(fid, bm.RbeamInf, 'int16');  % infinite curvature radius
    fwrite(fid, bm.pz, 'double');       % beam position z (m)
    fwrite(fid, bm.w0_pz, 'double');    % beam waist position z (m)
    fwrite(fid, bm.w0, 'double');       % beam waist radius (m)
    fwrite(fid, bm.zRay, 'double');     % Rayleigh distance (m)
    if     bm.PropType == 'In__to_Out'
      fwrite(fid, 'INSIDE__TO_OUTSIDE', 'char');
    elseif bm.PropType == 'Out_to_In_'
      fwrite(fid, 'OUTSIDE_TO_INSIDE_', 'char');
    else
      fwrite(fid, 'INSIDE__TO_INSIDE_', 'char');
    end
    fwrite(fid, bm.fr, 'double');       % focal ratio
    fwrite(fid, bm.diam, 'double');     % beam diameter (m)
    fclose(fid);

    if save_state_lam(1) == 0
      save_state_lam = bm.wl;
    else
      save_state_lam = [save_state_lam bm.wl];
    end
  end                   % if save_state ~= 0
end                     % function prop_savestate
