function [X_hat,Y_hat,y_hat] = StateCorrection(X,x_hat,param)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

R0       = param.earth.R0;
f        = param.earth.f;
g0       = param.earth.g0;
R0_INS_e = param.INS.Pos0_ECEF;

T = param.sim.T;

GPSTC_enabled = param.sensors.GPS.GPSTC_enabled;
ALT_enabled = param.sensors.ALT.ALT_enabled;

lat  = X(1);
long = X(2);
alt  = X(3);
V_c  = X(4:6);
q_bp = X(7:10);

% ---------
IMUbias = X(15:17);
IMUdrift = X(18:20);
% ---------

DR_hat         = x_hat(1:3);
DV_hat         = x_hat(4:6);
psi_hat        = x_hat(7:9);
nabla_hat      = x_hat(10:12);
epsilon_hat    = x_hat(13:15);

R_e_INS = LLA2ECEF(lat,long,alt);
D_ce_INS = DCM(2,-(lat+pi/2))*DCM(3,long);
R_e_hat = R_e_INS - D_ce_INS'*DR_hat;

LLA_hat = ECEF2LLA(R_e_hat);
lat_hat = LLA_hat(1);
long_hat = LLA_hat(2);
alt_hat = LLA_hat(3);

% Speed correction
V_hat = V_c - DV_hat;

% Oblate Earth
RN_hat = R0*(1-f*(2-3*(sin(lat_hat))^2)); %% or lat_GPS ??
RE_hat = R0*(1+f*(sin(lat_hat)^2));
Re = R0*(1-f*(sin(lat_hat))^2);
g_NED = g0*(1+0.0052884*(sin(lat_hat))^2)*(1-2*alt_hat/Re);

% Misalignment correction


% Total misalignment
dtheta_hat = [
    DR_hat(2)/(RE_hat+alt_hat)
    -DR_hat(1)/(RN_hat+alt_hat)
    -DR_hat(2)*tan(lat_hat)/(RE_hat+alt_hat)
    ];

phi_hat = psi_hat + dtheta_hat;

D_bp = quat2DCM(1,q_bp);
D_pt = eye(3)-skew(phi_hat);
D_bt_hat = D_bp*D_pt;
q_bt_hat = DCM2quat(D_bt_hat);
q_bt_hat = q_bt_hat/norm(q_bt_hat);

% -----------------------
% IMUbias_hat = IMUbias + nabla_hat;
% IMUdrift_hat = IMUdrift + epsilon_hat;
IMUbias_hat = IMUbias + nabla_hat;
IMUdrift_hat = IMUdrift + epsilon_hat;

X_NAV = [
    LLA_hat
    V_hat
    q_bt_hat
    0 %csi
    ];

X_IMU = [
    zeros(3,1)
    IMUbias_hat
    IMUdrift_hat
    ];

X_hat = [
    X_NAV
    X_IMU
    ];

% ----------------------
if GPSTC_enabled && false % <----------------------------------------------
    
    GPSbias = X(21);
    GPSbiasDot = X(22);
    
    biasGPS_hat    = x_hat(16);
    biasDotGPS_hat = x_hat(17);
    
    GPSbiasDot_hat = GPSbiasDot + biasDotGPS_hat;
    GPSbias_hat = GPSbias + biasGPS_hat;% + biasDotGPS_hat*T;
    X_GPS = [
        GPSbias_hat
        GPSbiasDot_hat
        ];
    X_hat = [
        X_hat
        X_GPS
        ];
end
% -----------------------
if ALT_enabled
    ALTbias = X(end);
    biasALT_hat = x_hat(end);
    ALTbias_hat = ALTbias + biasALT_hat;
    X_hat = [
        X_hat
        ALTbias_hat
        ];
end
% -----------------------

% Other outputs
D_ce_hat = DCM(2,-(lat_hat+pi/2))*DCM(3,long_hat);
Pos_NED  = D_ce_hat*(R_e_hat - R0_INS_e);
% R_INS_GPS_NED = D_ce_hat*R_e_hat;

euler_INS_GPS = DCM2euler(D_bt_hat,'ZYX');

% VN = V_hat(1);
% VE = V_hat(2);
% VD = V_hat(3);
% VG = sqrt(VN^2+VE^2);
% FPA = atan2(-VD,VG);
% TKA = atan2(VE,VN);

Y_hat = [
    R_e_INS
    Pos_NED
    euler_INS_GPS
    g_NED
%     VG
%     FPA
%     TKA
    ];

y_hat = [
    phi_hat
    psi_hat
    dtheta_hat
    ];
end