function [z_MAG,H_MAG,R_MAG,aux_MAG] = MAG_MM(MAG_data,X,U,param)

GPSTC_enabled = param.sensors.GPS.GPSTC_enabled;
ALT_enabled = param.sensors.ALT.ALT_enabled;

f  = param.earth.f;
R0 = param.earth.R0;

stdMAG = param.sensors.MAG.stdMAG;
stdMAGgain = param.sensors.MAG.stdMAGgain;

Bm_b = MAG_data(1:3);
% B_n  = MAG_data(4:6);

lat = X(1);
long = X(2);
alt = X(3);
q_bp = X(7:10);

% DRN = x(1);
% DRE = x(2);
% % DRD = x(3);
% psi = x(7:9);

RN = R0*(1-f*(2-3*(sin(lat))^2));
RE = R0*(1+f*(sin(lat)^2));

% theta = [DRE/(RE+alt) -DRN/(RN+alt) -tan(lat)*DRE/(RE+alt)]';
% D_ct = eye(3)-skew(theta);
% D_pc = eye(3)-skew(psi);
D_bp = quat2DCM(1,q_bp);
% D_bt = D_bp*D_pc*D_ct;

WMM_out = WMM(lat,long,alt);
Bc_n = WMM_out(1:3);

% z_MAG = Bm_b - D_bt*Bc_n;
z_MAG = Bm_b - D_bp*Bc_n;
% z_MAG = D_bp'*Bm_b - Bc_n;

Bmatrix = zeros(3);
Bmatrix(1,2) = 1/(RE+alt);
Bmatrix(2,1) = -1/(RN+alt);
Bmatrix(3,2) = -tan(lat)/(RE+alt);
H_MAG = D_bp*skew(Bc_n)*[Bmatrix zeros(3) eye(3) zeros(3,6)];
% H_MAG = skew(Bc_n)*[Bmatrix zeros(3) eye(3) zeros(3,6)];

if GPSTC_enabled
    H_MAG = [H_MAG zeros(3,2)];
end
if ALT_enabled
    H_MAG = [H_MAG zeros(3,1)];
end

R_MAG = eye(3)*(stdMAG*stdMAGgain)^2;

UB_MAG = 0;
aux_MAG = {'MAG',UB_MAG};

end