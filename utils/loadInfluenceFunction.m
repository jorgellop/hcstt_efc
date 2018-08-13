function influence_func = loadInfluenceFunction( filename, ac_spacing)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Read the influence function data from a FITS file
  info = fitsinfo(filename);
  infl  = fitsread(filename);

  [ldef, idef] = ismember('NAXIS1' , info.PrimaryData.Keywords(:, 1));
  ifnx = info.PrimaryData.Keywords{idef, 2};    % inf func number of pixels x

  [ldef, idef] = ismember('NAXIS2' , info.PrimaryData.Keywords(:, 1));
  ifny = info.PrimaryData.Keywords{idef, 2};    % inf func number of pixels y

  [ldef, idef] = ismember('P2PDX_M', info.PrimaryData.Keywords(:, 1));
  ifdx = info.PrimaryData.Keywords{idef, 2};    % inf func spacing x (m)

  [ldef, idef] = ismember('P2PDY_M', info.PrimaryData.Keywords(:, 1));
  ifdy = info.PrimaryData.Keywords{idef, 2};    % inf func spacing y (m)

  [ldef, idef] = ismember('C2CDX_M', info.PrimaryData.Keywords(:, 1));
  acdx = info.PrimaryData.Keywords{idef, 2};    % actuator spacing x (m)

  [ldef, idef] = ismember('C2CDY_M', info.PrimaryData.Keywords(:, 1));
  acdy = info.PrimaryData.Keywords{idef, 2};    % actuator spacing y (m)


ac_spacing0 = acdx/ifdx;
scale_factor = ac_spacing/ac_spacing0;
% ifnx_new = ifnx*scale_factor; 

  [cfx, cfy]   = meshgrid(-floor(ifnx/2):floor(ifnx/2),-floor(ifny/2):floor(ifny/2)); % meshgrids for influence function
  [cfx2, cfy2] = meshgrid(-floor(scale_factor*ifnx/4):floor(scale_factor*ifnx/4),-floor(scale_factor*ifny/4):floor(scale_factor*ifny/4)); % meshgrids for influence function


%[cfx_rs, cfy_rs] = meshgrid(-round(ifnx_new/2):round(ifnx_new/2),-round(ifnx_new/2):round(ifnx_new/2)); % meshgrids for influence function
influence_func = interp2(cfx*scale_factor,cfy*scale_factor,infl,cfx2,cfy2,'linear',0);
influence_func = influence_func/max(max(influence_func));

end

