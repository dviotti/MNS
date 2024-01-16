function AspNED = AspGT(U,param)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

R0        = param.earth.R0;
OmegaE    = param.earth.OmegaE;
f         = param.earth.f;
g0        = param.earth.g0;

lat      = U(1);
% long   = U(2);
h        = U(3);
VN       = U(4);
VE       = U(5);
VD       = U(6);
lat_dot  = U(7);
long_dot = U(8);
% h_dot =  U(9);
VN_dot   = U(10);
VE_dot   = U(11);
VD_dot   = U(12);

% Oblate Earth
% RE = R0*(1+f*(sin(lat))^2);
% RN = R0*(1-f*(2-3*(sin(lat))^2));
Re = R0*(1-f*(sin(lat))^2);
g_NED = g0*(1+0.0052884*(sin(lat))^2)*(1-2*h/Re);

% Position
% lat_dot = VN/(RN+h);
% long_dot = VE/((RE+h)*cos(lat));
% h_dot = -VD;

% Velocity
W1 = lat_dot;
W2 = (2*OmegaE+long_dot)*sin(lat);
W3 = (2*OmegaE+long_dot)*cos(lat);

AspN = VN_dot - W1*VD + W2*VE;
AspE = VE_dot - W2*VN - W3*VD;
AspD = VD_dot + W3*VE + W1*VN - g_NED;

AspNED = [
    AspN
    AspE
    AspD
    ];
end

