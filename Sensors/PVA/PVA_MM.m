function [z_PVA,H_PVA,R_PVA,aux] = PVA_MM(Z_PVA,X,param)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

PVAmodel = param.sensors.PVA.PVAmodel;

GPSTC_enabled = param.sensors.GPS.GPSTC_enabled;
ALT_enabled = param.sensors.ALT.ALT_enabled;

stdPVApos = param.sensors.PVA.stdPVApos;
stdPVAvel = param.sensors.PVA.stdPVAvel;
stdPVAatt = param.sensors.PVA.stdPVAatt;

Z_Pos = Z_PVA{1};
Z_Vel = Z_PVA{2};
Z_Att = Z_PVA{3};

lat  = X(1);
long = X(2);
alt  = X(3);
V_c  = X(4:6);
q_bp = X(7:10);

R_e = LLA2ECEF(lat,long,alt); % ECEF
D_ce = DCM(2,-(lat+pi/2))*DCM(3,long);
R_c = D_ce*R_e; % NEDt

z_PVApos = Z_Pos - R_c;
z_PVAvel = Z_Vel - V_c;

H_PVApos = [eye(3) zeros(3) zeros(3,9)];
H_PVAvel = [zeros(3) eye(3) zeros(3,9)];
H_PVAatt = [zeros(3,6) eye(3) zeros(3,6)];

% z_PVApos = H_PVApos*x_hat + randn(3,1)*stdPVApos;
% z_PVAvel = H_PVAvel*x_hat + randn(3,1)*stdPVAvel;
% z_PVAatt = H_PVAvel*x_hat + randn(3,1)*stdPVAatt;

R_PVApos = stdPVApos^2*eye(3);
R_PVAvel = stdPVAvel^2*eye(3);
R_PVAatt = stdPVAatt^2*eye(3);

switch PVAmodel
    case 1     
        H_PVA = H_PVApos;
        z_PVA = z_PVApos;
        R_PVA = R_PVApos;
    case 2
        H_PVA = H_PVAvel;
        z_PVA = z_PVAvel;
        R_PVA = R_PVAvel;
    case 3
        H_PVA = H_PVAatt;
        z_PVA = z_PVAatt;
        R_PVA = R_PVAatt;
    case 4
        H_PVA = [H_PVApos; H_PVAvel];
        z_PVA = [z_PVApos; z_PVAvel];
        R_PVA = blkdiag(R_PVApos,R_PVAvel);
    case 5
        H_PVA = [H_PVApos; H_PVAatt];
        z_PVA = [z_PVApos; z_PVAatt];
        R_PVA = blkdiag(R_PVApos,R_PVAatt);
    case 6
        H_PVA = [H_PVAvel; H_PVAatt];
        z_PVA = [z_PVAvel; z_PVAatt];
        R_PVA = blkdiag(R_PVAvel,R_PVAatt);
    case 7
        H_PVA = [H_PVApos; H_PVAvel; H_PVAatt];
        z_PVA = [z_PVApos; z_PVAvel; z_PVAatt];
        R_PVA = blkdiag(R_PVApos,R_PVAvel,R_PVAatt);
    otherwise
        error('Select a case')
end

if GPSTC_enabled
    H_PVA = [H_PVA zeros(size(H_PVA,1),2)];
end
if ALT_enabled
    H_PVA = [H_PVA zeros(size(H_PVA,1),1)];
end
UB = 0;
aux = {'PVA',UB};

end