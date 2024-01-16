function [z_GPS,H_GPS,R_GPS,aux_GPS] = GPSTC_MM(Z,X,U,param)

% GPStype = param.GPS.GPStype;
OmegaE = param.earth.OmegaE;
lightSpeed = param.earth.c;

% bias_reset    = GPS_reset(1);
% biasDot_reset = GPS_reset(2);
ALT_enabled = param.sensors.ALT.ALT_enabled;

% GPSnoise_rho = param.sensors.GPS.GPSnoise_rho;
GPSnoise_vel = param.sensors.GPS.GPSnoise_vel;
GPSarm = param.sensors.GPS.GPSarm;

stdPRgain = param.sensors.GPS.stdPRgain;
stdDPRgain = param.sensors.GPS.stdDPRgain;
DPRon = param.sensors.GPS.DPRon;

atm_corr_enabled = param.sensors.GPS.atm_corr_enabled;
Erot_corr_enabled = param.sensors.GPS.Erot_corr_enabled;
alpha = param.sensors.GPS.ion_corr.alpha;
beta = param.sensors.GPS.ion_corr.beta;

d = param.sensors.GPS.noise.d;
B = param.sensors.GPS.noise.B;

% L = param.sim.L;

lat  = X(1);
long = X(2);
alt  = X(3);
V_c  = X(4:6);
q_bp = X(7:10);

% psi_hat = x_hat(7:9);
% bias = x_hat(16);
% bias = X(21); <----------------------------------------------------------
% biasDot = X(22); <-------------------------------------------------------

omega = U(4:6);

R_INS_e = LLA2ECEF(lat,long,alt);
D_ce_INS = DCM(2,-(lat+pi/2))*DCM(3,long);

nZ = numel(Z);
switch nZ
    case 5
        SVID = Z{1};
        R_SV_eT = Z{2};
        rho_m = Z{3};
        URA = Z{4};
        CN0 = Z{5};
        GPSvel_enabled = false;
    case 7
        SVID = Z{1};
        R_SV_eT = Z{2};
        rho_m = Z{3};
        URA = Z{4};
        CN0 = Z{5};
        V_SV_eT = Z{6};
        vel_m = Z{7};
        GPSvel_enabled = true;
    otherwise
        error('Z_GPSTC has worng number of elements')
end

[p,q] = size(rho_m);
if p==1
    rho_m = rho_m';
elseif ~(p>=1 && q==1)
    error('Peseudo range must be a mx1 vector')
end

visibleSVs = numel(rho_m);

D_bp = quat2DCM(1,q_bp);
% D_pc = eye(3) - skew(psi_hat);
% D_bc = D_bp*D_pc;
% LOS_b = D_bc*LOS_c;

R_Rx_e = R_INS_e + (D_bp*D_ce_INS)'*GPSarm;
[LOS_eT,~] = GPSlinear(R_SV_eT,R_Rx_e);
LOS_c = D_ce_INS*LOS_eT';

D = 60;
tow = 387000;
rho_corr = NaN(visibleSVs,1);
R_SV_eR = NaN(3,visibleSVs);
var = NaN(visibleSVs,1);
for i=1:visibleSVs
    
    az_c = wrapTo2Pi(atan2(LOS_c(2,i),LOS_c(1,i)));
    hip_c = sqrt(LOS_c(1,i)^2+LOS_c(2,i)^2);
    el_c = wrapToPi(atan2(-LOS_c(3,i),hip_c));
    
    % Atmosphere corrections
    if atm_corr_enabled
        [~, dIon,~,varIono] = IonosphericDelay(tow,lat,long,az_c,el_c,alpha,beta);
        [dTr,stdTr] = TropoCorr(lat,alt,el_c,D,param);
        rho_corr(i) = rho_m(i) - (dIon + dTr);% - bias;
        % Atm variances
        stdNoise = sqrt(d*B/(2*CN0(i)))*293;
        stdMP = 0.13 + 0.53*exp(-el_c*180/pi/10);
        varAir = stdMP^2 + stdNoise^2;
        % varAir = 25; % DO-316 J.3
        varAtm = varIono + stdTr^2 + varAir;
    else
        rho_corr(i) = rho_m(i);% - bias; <---------------------------------
        varAtm = 0;
    end
    
    % Earth rotation correction
    if Erot_corr_enabled
        angle = OmegaE*rho_corr(i)/lightSpeed; % should be angle = OmegaE*(rho_m(i)-bias)/lightSpeed;
        R_SV_eR(:,i) = DCM(3,angle)*R_SV_eT(:,i);
    else
        R_SV_eR(:,i) = R_SV_eT(:,i);
    end
    
    % Pseudorange variance
    var(i) = URA(i) + varAtm;
end

% Pseudo-distance
[LOS_eR,rho_c] = GPSlinear(R_SV_eR,R_Rx_e);
H_rho_pos = LOS_eR*D_ce_INS';
H_rho_att = H_rho_pos*skew(D_bp'*GPSarm);%zeros(visibleSVs,3);
H_rho = [H_rho_pos zeros(visibleSVs,3) H_rho_att zeros(visibleSVs,6) ones(visibleSVs,1) zeros(visibleSVs,1)];
% z_rho = rho_m - rho_c;% - bias_reset*ones(visibleSVs,1);
z_rho = rho_corr - rho_c;
% R_rho = eye(visibleSVs)*GPSnoise_rho^2;
R_rho = diag(var.*stdPRgain^2);
% Velocity
if GPSvel_enabled && DPRon
    H_rhoDot_vel = LOS_eR*D_ce_INS';
    H_rhoDot_psi = H_rhoDot_vel*skew(omega)*skew(D_bp'*GPSarm);
    H_rhoDot_drift = H_rhoDot_vel*skew(D_bp'*GPSarm);
    H_rhoDot = [zeros(visibleSVs,3) H_rhoDot_vel H_rhoDot_psi zeros(visibleSVs,3) ...
        H_rhoDot_drift zeros(visibleSVs,1) ones(visibleSVs,1)];
    Mvel_INS = LOS_eR*(V_SV_eT - D_ce_INS'*V_c);
    vel_INS = diag(Mvel_INS);
    z_vel = vel_m - vel_INS;% - biasDot*ones(visibleSVs,1); <--------------
    R_vel = eye(visibleSVs)*(GPSnoise_vel*stdDPRgain)^2;
    % GPS observation model
    H_GPS = [H_rho; H_rhoDot];
    z_GPS = [z_rho; z_vel];
    R_GPS = blkdiag(R_rho,R_vel);
else
    H_GPS = H_rho;
    z_GPS = z_rho;
    R_GPS = R_rho;
end
% H_GPS = H_GPS*L;

if ALT_enabled
    if GPSvel_enabled && DPRon
        H_GPS = [H_GPS zeros(2*visibleSVs,1)];
    else
        H_GPS = [H_GPS zeros(visibleSVs,1)];
    end
end

% dof_z_GPS = min(4,visibleSVs);
% if GPSvel_enabled
%     dof_z_GPS = 2*dof_z_GPS;
% end
% % dof_z_GPS = numel(rho_m);

UB_GPS = zeros(numel(z_GPS),1);%sum((1./rho_corr).^2);

aux_GPS = {'GPSTC',UB_GPS,SVID};

end