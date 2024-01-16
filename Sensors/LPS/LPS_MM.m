function [z_LPS,H_LPS,R_LPS,aux] = LPS_MM(Z,X,U,param)

stdLPS = param.sensors.LPS.stdLPS;
Ra_e = param.sensors.LPS.Ra_e;
stdLPSgain = param.sensors.LPS.stdLPSgain;

GPSTC_enabled = param.sensors.GPS.GPSTC_enabled;
ALT_enabled = param.sensors.ALT.ALT_enabled;

dist_m = Z{1};
tracking = Z{2};

lat = X(1);
long = X(2);
alt = X(3);

dist_tracking = dist_m(tracking);
Ra_tracking = Ra_e(:,tracking);

n = numel(dist_tracking);

R_rec = LLA2ECEF(lat,long,alt);
[LOS_e,dist_c] = GPSlinear(Ra_tracking,R_rec);

z_LPS = dist_tracking' - dist_c;

D_ne = DCM(2,-(lat+pi/2))*DCM(3,long);

H_LPS = [LOS_e*D_ne' zeros(n,12)];
if GPSTC_enabled
    H_LPS = [H_LPS zeros(n,2)];
end
if ALT_enabled
    H_LPS = [H_LPS zeros(n,1)];
end

R_LPS = eye(n)*(stdLPS*stdLPSgain)^2;

UB = 1./(dist_c.^2);
aux = {'LPS',UB};

end