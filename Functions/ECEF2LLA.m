function LLA = ECEF2LLA(R)

% Parameters
% inv_a         = 1.56785594288739799723e-7;
inv_a2        = 2.45817225764733181057e-14;
l             = 3.34718999507065852867e-3;
l2            = 1.12036808631011150655e-5;
one_e2        = 9.93305620009858682943e-1;
one_e2_div_a2 = 2.44171631847341700642e-14;
% one_e2_div_b  = 1.56259921876129741211e-7;
Hmin          = 2.25010182030430273673e-14;
inv_cbcrt2    = 7.93700525984099737380e-1;
% a2_div_c      = 7.79540464078689228919e7;
% b2_div_c2     = 1.48379031586596594555e2;

x = R(1);
y = R(2);
z = R(3);

w2 = x^2+y^2;

m = w2*inv_a2;
n = z^2*one_e2_div_a2;
p = (m+n-4*l2)/6;

G = m*n*l2;
H = 2*p^3+G;

if H<Hmin
    error('H < Hmin')
end

aux = H+G+2*sqrt(H*G);
C = nthroot(aux,3)*inv_cbcrt2;

i = -(2*l2+m+n)/2;
P=p^2;
beta = i/3-C-P/C;
k = l2*(l2-m-n);

t = sqrt(sqrt(beta^2-k)-(beta+i)/2) - sign(m-n)*sqrt(abs(beta-i)/2);

F = t^4+2*i*t^2+2*l*(m-n)*t+k;
dFdt = 4*t^3*4*i*t+2*l*(m-n);
Dt = -F/dFdt;

u = t+Dt+l;
v = t+Dt-l;

w = sqrt(w2);

lat = atan2(z*u,w*v);

Dw = w*(1-1/u);
Dz = z*(1-one_e2/v);

h = sign(u-1)*sqrt(Dw^2+Dz^2);

long = atan2(y,x);

LLA = [
    lat
    long
    h
    ];

end





