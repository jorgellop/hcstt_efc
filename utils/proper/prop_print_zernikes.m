function prop_print_zernikes(numz)
%        prop_print_zernikes(numz)
% Print to the screen the first numz of Noll-ordered Zernike polynomials
% for an unobscured circular aperture.
% numz = number of Zernike polynomials to print (1 to numz)

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Revision history:
% 2005 Feb     jek
% 2016 May 02  gmg  Translated to Matlab

  list = prop_noll_zernikes(numz);
  for iz = 1 : numz
    fprintf(1, '%8d  =  %s\n', iz, char(list(iz)));
  end
end                     % prop_print_zernikes
