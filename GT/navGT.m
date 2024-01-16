function [Xdot,Y] = navGT(t,X,U,param)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

R0        = param.earth.R0;
% OmegaE    = param.earth.OmegaE;
f         = param.earth.f;
% g0        = param.earth.g0;
Pos0_ECEF = param.navGT.Pos0_ECEF;

lat  = X(1);
long = X(2);
h    = X(3);
VN   = X(4);
VE   = X(5);
VD   = X(6);

VN_dot = U(1);
VE_dot = U(2);
VD_dot = U(3);

% Oblate Earth
RE = R0*(1+f*(sin(lat))^2);
RN = R0*(1-f*(2-3*(sin(lat))^2));
% Re = R0*(1-f*(sin(lat))^2);
% g_NED = g0*(1+0.0052884*(sin(lat))^2)*(1-2*h/Re);

lat_dot = VN/(RN+h);
long_dot = VE/((RE+h)*cos(lat));
h_dot = -VD;

Xdot = [
    lat_dot
    long_dot
    h_dot
    VN_dot
    VE_dot
    VD_dot
    ];

Pos_ECEF = LLA2ECEF(lat,long,h);
xyz_ECEF = Pos_ECEF-Pos0_ECEF;
D_ve = DCM(2,-(lat+pi/2))*DCM(3,long);
xyz_NED  = D_ve*xyz_ECEF;

Y = [
    lat_dot
    long_dot
    h_dot
    VN_dot
    VE_dot
    VD_dot
    Pos_ECEF
    xyz_NED
    ];

end

