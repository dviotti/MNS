%% ======================== INS INITIALIZATION ============================

% LLA0
if GPSTC_enabled
    R_SV_e = Z_GPSTC{2};
    rho_m = Z_GPSTC{3};
    URA_m = Z_GPSTC{4};
    if Erot_corr_enabled
        [~,R0_rec,b0] = GPSsolution(R_SV_e,rho_m,URA_m,param);
    else
        [R0_rec,b0] = GPS_LS(R_SV_e,rho_m,URA_m);
    end
    GPS_LLA = ECEF2LLA(R0_rec);
    lat0 = GPS_LLA(1);
    long0 = GPS_LLA(2);
    alt0 = GPS_LLA(3);% + paramAC.GPSarm(3);
    % TODO: LS solution
    V0_N = 0;
    V0_E = 0;
    V0_D = 0;
elseif GPSLC_enabled
    GPS_LLA = Z_GPSLC{1};
    lat0 = GPS_LLA(1);
    long0 = GPS_LLA(2);
    alt0 = GPS_LLA(3);% + paramAC.GPSarm(3);
    if numel(Z_GPSLC)==2
        GPS_vel = Z_GPSLC{2};
        V0_N = GPS_vel(1);
        V0_E = GPS_vel(2);
        V0_D = GPS_vel(3);
    else
        V0_N = 0;
        V0_E = 0;
        V0_D = 0;
    end
else
    lat0 = lat_GT(1);
    long0 = long_GT(2);
    alt0 = alt_GT(3);
    V0_N = V_l_GT(1,1);
    V0_E = V_l_GT(2,1);
    V0_D = V_l_GT(3,1);
end

% Velocity related signals
%     V0_G = sqrt(V0_N^2+V0_E^2); %Ground speed
%     FPA0 = atan2(-V0_D,V0_G);
%     TKA0 = atan2(V0_E,V0_N);

% Origin in ECEF and computed NED
R0_INS_e = LLA2ECEF(lat0,long0,alt0);
D0_ce = DCM(2,-(lat0+pi/2))*DCM(3,long0);
R0_INS_c = D0_ce*R0_INS_e;

% Oblate Earth
RE = R0*(1+f*(sin(lat0))^2);
RN = R0*(1-f*(2-3*(sin(lat0))^2));
Re = R0*(1-f*(sin(lat0))^2);
g_NED0 = g0*(1+0.0052884*(sin(lat0))^2)*(1-2*alt0/Re);

% ------------------- INITIAL ALIGNMENT: TRIAD ------------------------
t0 = 30;
if TRIAD_on
    % Median of 30 sec data
    k_INS = round(t0/T)+1;
    k_IMU = round(t0/T_IMU)+1;
    k_MAG = round(t0/param.sensors.MAG.T_MAG)+1;
    if strcmp(IMU_DATAtype,'SIM')
        Asp0_b = median(dAsp_b_m(:,1:k_IMU),2)/T_IMU;
        omega0_bi_b = median(dOmega_bi_b_m(:,1:k_IMU),2)/T_IMU;
        B_bm = median(MAG_data(:,1:k_MAG),2);
    elseif strcmp(IMU_DATAtype,'FT')
        Asp0_b = median(Asp(:,Dt),2)/T;
        omega0_bi_b = median(Omega(:,Dt),2)/T;
    end
    
    g0_b = -Asp0_b;
    g0 = norm(g0_b); % Initial gravity magnitude
    
    Re = R0*(1-f*(sin(lat0))^2);
    g_l = g0*(1+0.0052884*(sin(lat0))^2)*(1-2*alt0/Re);
    Asp0_l = [0 0 -g_l]';
    
    if TRIAD_MAG
        V2_b = B_bm;
        B_c = WMM(lat0,long0,alt0);
        V2_l = B_c(1:3);
    else
        V2_b = omega0_bi_b;
        OmegaE_e = [0 0 OmegaE]';
        OmegaE_l = D0_ce*OmegaE_e;
        V2_l = OmegaE_l;
    end
    
    % ### TRIAD ###
    % Body Frame
    v1_b = Asp0_b;
    v1_b = v1_b/norm(v1_b);
    v2_b = skew(Asp0_b)*V2_b;
    v2_b = v2_b/norm(v2_b);
    v3_b = skew(Asp0_b)*skew(Asp0_b)*V2_b;
    v3_b = v3_b/norm(v3_b);
    Y_b = [v1_b v2_b v3_b];
    % NED
    v1_l = Asp0_l;
    v1_l = v1_l/norm(v1_l);
    v2_l = skew(Asp0_l)*V2_l;
    v2_l = v2_l/norm(v2_l);
    v3_l = skew(Asp0_l)*skew(Asp0_l)*V2_l;
    v3_l = v3_l/norm(v3_l);
    X_l = [v1_l v2_l v3_l];
    
    [q0_bl,euler0,D0_bv] = TRIAD(X_l,Y_b);
    
else
    q0_bl = q_bl_GT(:,1);
    D0_bv = quat2DCM(1,q0_bl);
    euler0 = DCM2euler(D0_bv,'ZYX');
end

% ---------------------------------------------------------------------

X0 = [
    lat0
    long0
    alt0
    V0_N
    V0_E
    V0_D
    q0_bl
    0
    zeros(3,1)
    zeros(6,1)
    ];
%     if GPSTC_enabled <---------------------------------------------------
%         X0_GPS = [b0; 100];%112
%         X0 = [X0; X0_GPS];
%     end
if ALT_enabled
    SF = (Z_ALT(1) - alt0)/alt0;
    X0 = [X0; SF];
    %         bias = (Z_ALT(1) - alt0);
    %         X0 = [X0; bias];
end
%     Y0 = [
%         Y_GT(1,1:9)'
%         g_NED0
% %         V0_G
% %         FPA0
% %         TKA0
%         zeros(6,1)
%         ];
Y0 = [
    R0_INS_e
    R0_INS_c
    euler0
    g_NED0
    %         V0_G
    %         FPA0
    %         TKA0
    zeros(6,1)
    ];

param.INS = struct('K',Kvert,'Pos0_ECEF',R0_INS_e);

X(:,1) = X0;
Y(:,1) = Y0;
%     V_G = V0_G;
