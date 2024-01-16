function R = LLA2ECEF(lat,long,h)

% Parameters
one_e2        = 9.93305620009858682943e-1;
a2_div_c      = 7.79540464078689228919e7;
b2_div_c2     = 1.48379031586596594555e2;

N = a2_div_c/sqrt(cos(lat)^2+b2_div_c2);
d = (N+h)*cos(lat);

x = d*cos(long);
y = d*sin(long);
z = (N*one_e2+h)*sin(lat);

R = [x y z]';

end