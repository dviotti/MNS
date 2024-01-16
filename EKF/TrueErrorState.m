function [x,y] = TrueErrorState(X,X_GT,Y_GT,param)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

R0 = param.earth.R0;
f  = param.earth.f;

IMU_DATAtype = param.sensors.IMU.IMU_DATAtype;

lat  = X(1);
long = X(2);
alt  = X(3);
V_c  = X(4:6);
q_bp = X(7:10);

lat_GT  = X_GT(1);
long_GT = X_GT(2);
alt_GT  = X_GT(3);
V_l_GT  = X_GT(4:6);
q_bt_GT = X_GT(7:10);


if strcmp(IMU_DATAtype,'FT')
    GPSarm = param.AC.GPSarm;
    % POSRx_e_GT = Y_GT(1:3);
    POS_INS_e_GT = Y_GT(1:3);
    euler = Y_GT(7:9);
    % ----
    D_bl_GT = euler2DCM(euler,'ZYX');
    D_le_GT = DCM(2,-(lat_GT+pi/2))*DCM(3,long_GT);
    R_e_INS = LLA2ECEF(lat,long,alt); % ECEF
    % DR_l_GT = D_le_GT*(R_e_INS - POSRx_e_GT) + D_bl_GT'*GPSarm;
    DR_l_GT = D_le_GT*(R_e_INS - POS_INS_e_GT) + D_bl_GT'*GPSarm;
    % ----
else
    % True attitude
    D_bl_GT = quat2DCM(1,q_bt_GT);
    
    % DR
    R_e_GT = LLA2ECEF(lat_GT,long_GT,alt_GT); % ECEF
    D_le_GT = DCM(2,-(lat_GT+pi/2))*DCM(3,long_GT);
    R_e_INS = LLA2ECEF(lat,long,alt); % ECEF
%     D_ce = DCM(2,-(lat+pi/2))*DCM(3,long);
    DR_l_GT = D_le_GT*(R_e_INS - R_e_GT); % NEDt
%     DR_l_GT = D_ce*(R_e_INS - R_e_GT); % NEDt
%     DR_l_GT = D_ce*R_e_INS - D_le_GT*R_e_GT; % NEDt
end

RN_GT = R0*(1-f*(2-3*(sin(lat_GT))^2));
RE_GT = R0*(1+f*(sin(lat_GT)^2));

dtheta_GT = [
    DR_l_GT(2)/(RE_GT+alt_GT)
    -DR_l_GT(1)/(RN_GT+alt_GT)
    -DR_l_GT(2)*tan(lat_GT)/(RE_GT+alt_GT)
    ];

% DV
% D_ct = eye(3)-skew(dtheta_GT);
% DV_l_GT = D_ct'*V_c - V_l_GT;
DV_l_GT = V_c - V_l_GT;

% psi
D_bt = D_bl_GT;
D_bp = quat2DCM(1,q_bp);
D_pt = D_bp'*D_bt;
phi_GT = skew(eye(3)-D_pt);

psi_GT = phi_GT - dtheta_GT;

% -----
x = [
    DR_l_GT
    DV_l_GT
    psi_GT
    ];

y = [
    phi_GT
    psi_GT
    dtheta_GT
    ];

end