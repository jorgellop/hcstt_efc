function zernikeMat = readZernikeMatrix()
% line1 = ['Optic','  ','tip','  ','tilt','  ','tilt','  ','defoc','  ','astig','  ','astig','  ','coma','  ','coma','  ','trefo','  ','trefo','  ','spher1','  ','trefo2','  ','spher2','  ','Z_6^6','  ','Z_7^3','  ','Z_7^7','  ','Z_8^4'];
% txtfile = fopen('LowOrderAberrations.txt','w');
% format = '%110s \r';% %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s';
% fprintf(txtfile,format,line1);

% mat=[0,0,0;...
%     0,1,1;
%     0,0,1];
% format = '%2u %2u %2u \r\n';
% fprintf(txtfile,format,mat');

fid = fopen('utils\LowOrderAberrations_test.txt', 'rt');
C = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1, 'CollectOutput', true);
fclose(fid);
optnames = C{1};
mat0 = C{2};
% data_name3 = C{3};
% zernikeMat_presc2DM1 = mat0(1:13)
% zernikeMat_DM12Image = mat0(14:end)

mat0sum = sum(mat0');
opt = find(mat0sum~=0);
optpad = padarray(opt,[0,17-length(opt)],'post');
zernikeMat = zeros(1+length(opt),17);
zernikeMat(1,:) = optpad;
zernikeMat(2:end,:) = mat0(opt,:);

end