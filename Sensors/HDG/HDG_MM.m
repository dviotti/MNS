function [z_HDG,H_HDG,R_HDG,aux_HDG] = HDG_MM(HDG_data,X,U,param)

GPSTC_enabled = param.sensors.GPS.GPSTC_enabled;
ALT_enabled = param.sensors.ALT.ALT_enabled;

stdHDG = param.sensors.HDG.stdHDG;
stdHDGgain = param.sensors.HDG.stdHDGgain;

HDGmag_m = HDG_data;

lat = X(1);
long = X(2);
alt = X(3);
q_bc = X(7:10);

D_bc = quat2DCM(1,q_bc);
% HDG_c = atan2(D_bc(3,2),D_bc(3,1));
HDGtrue_c = atan2(D_bc(1,2),D_bc(1,1));

MagVar = WMM(lat,long,alt);
MagDecl = MagVar(7);
HDGmag_c = HDGtrue_c - MagDecl; %<- check sign

z_HDG = HDGmag_m - HDGmag_c;

H_HDG = [zeros(1,8) 1 zeros(1,6)];

if GPSTC_enabled
    H_HDG = [H_HDG zeros(1,2)];
end
if ALT_enabled
    H_HDG = [H_HDG 0];
end

R_HDG = (stdHDGgain*stdHDG)^2;

UB = 0;
aux_HDG = {'HDG',UB};

end