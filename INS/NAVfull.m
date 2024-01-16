function [Xdot,Y] = NAVfull(t,X,U,param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

R0          = param.earth.R0;
OmegaE      = param.earth.OmegaE;
f           = param.earth.f;
g0          = param.earth.g0;
Pos0_ECEF   = param.INS.Pos0_ECEF;
% Pos0_LLA    = param.sim.Pos0_LLA;
K           = param.INS.K;
IMUarm      = param.sensors.IMU.IMUarm;
wn          = param.sensors.IMU.wn;

lat = X(1);
long = X(2);
h = X(3);
VN = X(4);
VE = X(5);
VD = X(6);
q_bv = X(7:10);
csi = X(11);
arm = X(12:14);

% Oblate Earth
RE = R0*(1+f*(sin(lat))^2);
RN = R0*(1-f*(2-3*(sin(lat))^2));
Re = R0*(1-f*(sin(lat))^2);
g_NED = g0*(1+0.0052884*(sin(lat))^2)*(1-2*h/Re);

% AspN_b = U(1);
% AspE_b = U(2);
% AspD_b = U(3);
Asp_b_IMU = U(1:3);%*g_NED;
omega_bi_b = U(4:6);
altimiter = U(7);

% Arm compensation
arm_dot = wn*(omega_bi_b - arm);
alpha = omega_bi_b - arm;
Asp_b_CG = Asp_b_IMU - (skew(alpha)*IMUarm + skew(omega_bi_b)*(skew(omega_bi_b)*IMUarm));

q_bv = q_bv/norm(q_bv);
D_bv = quat2DCM(1,q_bv);
D_vb = D_bv';
Asp_v = D_vb*Asp_b_CG;

AspN = Asp_v(1);
AspE = Asp_v(2);
AspD = Asp_v(3);

% Vertical channel stabilization 
k1 = K(1);
k2 = K(2);
k3 = K(3);
% h0 = Pos0_LLA(3);
err = h - altimiter;
csi_dot = k3*err;

% Position
lat_dot = VN/(RN+h);
long_dot = VE/((RE+h)*cos(lat));
h_dot = -VD - k1*err;

% Velocity
W1 = lat_dot;
W2 = (2*OmegaE+long_dot)*sin(lat);
W3 = (2*OmegaE+long_dot)*cos(lat);

VN_dot = AspN + W1*VD - W2*VE;
VE_dot = AspE + W2*VN + W3*VD;
VD_dot = AspD - W3*VE - W1*VN + g_NED + k2*err + csi;

Omega_b = [
    0 -omega_bi_b'
    omega_bi_b -skew(omega_bi_b)
    ];

omega_vi_v = [
    (OmegaE+long_dot)*cos(lat)
    -lat_dot
    -(OmegaE+long_dot)*sin(lat)
    ];
Omega_v = [
    0 omega_vi_v'
    -omega_vi_v -skew(omega_vi_v)
    ];

q_bv_dot = 1/2*(Omega_b+Omega_v)*q_bv;

Xdot = [
    lat_dot
    long_dot
    h_dot
    VN_dot
    VE_dot
    VD_dot
    q_bv_dot
    csi_dot
    arm_dot
    ];

Pos_ECEF = LLA2ECEF(lat,long,h);
D_ve = DCM(2,-(lat+pi/2))*DCM(3,long);
Pos_NED  = D_ve*(Pos_ECEF-Pos0_ECEF);

euler = DCM2euler(D_bv,'ZYX');

% VG = sqrt(VN^2+VE^2);
% FPA = atan2(-VD,VG);
% TKA = atan2(VE,VN);

omega_bv_b = omega_bi_b - D_bv*omega_vi_v;

Y = [
    Pos_ECEF
    Pos_NED
    euler
    g_NED
    Asp_b_CG
    omega_bv_b
%     VG
%     FPA
%     TKA
%     Asp_b_CG
%     Asp_b_IMU
    ];

end

