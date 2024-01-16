function [R_SVvis_e,V_SVvis_E_e,range_vis,nu_vis,az_b,el_b,az_l,el_l,isVisible] ...
    = SVvisible(t,X,U,R_SV_eR,V_SV_E_eR,param)

c = param.earth.c;
f_L1 = param.SV.f_L1;
el_mask = param.SV.elevation_rad;
SV_dt = param.SV.SV_dt;
N_SV = param.SV.N_SV;
idx0 = param.SV.idx0;
OmegaE = param.earth.OmegaE;

a_Rx_b = param.SV.a_Rx_b;

lambda_L1 = c/f_L1;

% eR: ECEF coordinate system at the time of receprion of the SV signal
% eT: ECEF coordinate system at the time of transmission of the SV signal

lat = X(1);
long = X(2);
h = X(3);
V_rec_E_l = X(4:6);
q_rec_l2b = X(7:10);

omega_BI_b = U(4:6);

% S_E to S_NED -> D_ve
D_leR = DCM(2,-(lat+pi/2))*DCM(3,long);
D_bl = quat2DCM(1,q_rec_l2b);
D_beR = D_bl*D_leR;

% R_rec_eR
R_eR = LLA2ECEF(lat,long,h);
a_Rx_eR = D_beR'*a_Rx_b;
R_rec_eR = R_eR + a_Rx_eR;

% V_rec_E_e
OmegaE_l = [OmegaE*cos(lat) 0 -OmegaE*sin(lat)]';
omega_BE_l =  D_bl'*omega_BI_b - OmegaE_l;
V_rec_E_eR = D_leR'*(V_rec_E_l + skew(omega_BE_l)*(D_bl'*a_Rx_b));

az_l = NaN(1,N_SV); % Azimute wrt local frame
el_l = NaN(1,N_SV); % Elevation wrt local frame
az_b = NaN(1,N_SV); % Azimute wrt body frame
el_b = NaN(1,N_SV); % Elevation wrt body frame

% mag_V_rel_LOS = zeros(1,N_SV);
range = zeros(1,N_SV); % range
nu = zeros(1,N_SV); % Doppler effect

R_SVj_eT = zeros(3,N_SV);
V_SVj_E_eT = zeros(3,N_SV);
isVisible = zeros(1,N_SV);
for j=1:N_SV
    
    i_corr = round(t/SV_dt)+idx0;
    R_SVj_eT(:,j) = R_SV_eR(i_corr,:,j)';
    
    Dt = 0;
    
    R_rec2SVj_eR = DCM(3,OmegaE*Dt)*R_SVj_eT(:,j) - R_rec_eR;
    range(j) = norm(R_rec2SVj_eR);
    
    Dt_prev = Dt;
    Dt = range(j)/c;
    i = 0;
    while abs(Dt-Dt_prev)>1e-6
        
        t_corr = t - Dt;
        i_corr = round(t_corr/SV_dt)+idx0;
        R_SVj_eT(:,j) = R_SV_eR(i_corr,:,j)';
        
        R_rec2SVj_eR = DCM(3,OmegaE*Dt)*R_SVj_eT(:,j) - R_rec_eR;
        range(j) = norm(R_rec2SVj_eR);
        
        Dt_prev = Dt;
        Dt = range(j)/c;
        
        if i>10
            warning('SV solution not found')
            break
        end
        i = i + 1;
    end
    
    LOS_eR = R_rec2SVj_eR/range(j);
    LOS_l = D_leR*LOS_eR;

    az_l(j) = wrapTo2Pi(atan2(LOS_l(2),LOS_l(1)));
    hip_l = sqrt(LOS_l(1)^2+LOS_l(2)^2);
    el_l(j) = wrapToPi(atan2(-LOS_l(3),hip_l));

    LOS_b = D_bl*LOS_l;

    az_b(j) = wrapTo2Pi(atan2(LOS_b(2),LOS_b(1)));
    hip_b = sqrt(LOS_b(1)^2+LOS_b(2)^2);
    el_b(j) = wrapToPi(atan2(-LOS_b(3),hip_b));

    if el_l(j)>el_mask && el_b(j)>0
        isVisible(j) = 1;
    else
        continue
    end
    
    V_SVj_E_eT(:,j) = V_SV_E_eR(i_corr,:,j)';
    V_rel_E_e = DCM(3,OmegaE*Dt)*V_SVj_E_eT(:,j) - V_rec_E_eR;
    mag_V_rel_LOS_eR = V_rel_E_e'*LOS_eR;
    %  Vj_rel_LOS = mag_V_rel_LOS*LOS_e;
    
    nu(j) = -mag_V_rel_LOS_eR/lambda_L1;
    
    % V_rel_LOS_ort = V_rel_e - Vj_rel_LOS;
    % w_LOS2e_e = V_rel_LOS_ort/norm(dist);
    %
    % D_lv = DCM(2,gamma(j))*DCM(3,chi(j));
    % D_le = D_lv*D_ve;
    % w_LOS2NED_l = D_le*w_LOS2e_e;
    %
    % gamma_dot(j) = w_LOS2NED_l(2);
    % chi_dot(j) = w_LOS2NED_l(3)/cos(gamma(j));

end

% visibleSVs = sum(isVisible);
R_SVvis_e = R_SVj_eT(:,isVisible==1);
V_SVvis_E_e = V_SVj_E_eT(:,isVisible==1);
% chi_vis = chi(isVisible==1);
% gamma_vis = gamma(isVisible==1);
range_vis = range(isVisible==1);
nu_vis = nu(isVisible==1);
% vel_vis = -nu_vis*lambda_L1;

end

