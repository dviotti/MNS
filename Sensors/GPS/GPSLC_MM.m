function [z_GPS,H_GPS,R_GPS,dof_z_GPS,UB] = GPSLC_MM(Z,X,x_hat,param)

stdGPSpos = param.sensors.GPSLC.stdGPSLCpos;
stdGPSvel = param.sensors.GPSLC.stdGPSLCvel;
stdGPSLCgain = param.sensors.GPSLC.stdGPSLCgain;

ALT_enabled = param.sensors.ALT.ALT_enabled;

lat  = X(1);
long = X(2);
alt  = X(3);
V_c  = X(4:6);
        
R_INS_e = LLA2ECEF(lat,long,alt);
% D_ce_INS = DCM(2,-(lat+pi/2))*DCM(3,long);

nZ = numel(Z);
switch nZ
    case 1
        GPS_LLA = Z{1};
        GPSvel_enabled = false;
    case 2
        GPS_LLA = Z{1};
        GPS_Vel = Z{2};
        GPSvel_enabled = true;
    otherwise
        error('Z_GPSLC has worng number of elements')
end

GPS_lat = GPS_LLA(1);
GPS_long = GPS_LLA(2);
GPS_alt = GPS_LLA(3);

Rrec_e = LLA2ECEF(GPS_lat,GPS_long,GPS_alt);
D_ce_GPS = DCM(2,-(GPS_lat+pi/2))*DCM(3,GPS_long);
z_R_GPS_c = D_ce_GPS*(R_INS_e - Rrec_e);
z_GPS = z_R_GPS_c;
H_GPS = [eye(3) zeros(3,12)];

R_GPSpos(3) = (2*stdGPSpos*stdGPSLCgain)^2;
R_GPSpos(2) = (stdGPSpos*stdGPSLCgain)^2;
R_GPSpos(1) = (stdGPSpos*stdGPSLCgain)^2;
R_GPSpos = diag(R_GPSpos);
R_GPS = R_GPSpos;

dof_z_GPS = 3;

if GPSvel_enabled
    z_V_GPS_c = V_c - GPS_Vel;
    z_GPS = [z_R_GPS_c; z_V_GPS_c];
    H_GPS = [eye(6) zeros(6,9)];
    
    R_GPSvel(3) = (2*stdGPSvel*stdGPSLCgain)^2;
    R_GPSvel(2) = (stdGPSvel*stdGPSLCgain)^2;
    R_GPSvel(1) = (stdGPSvel*stdGPSLCgain)^2;
    R_GPSvel = diag(R_GPSvel);
    R_GPS = blkdiag(R_GPS,R_GPSvel);
    
    dof_z_GPS = 6;
end

if ALT_enabled
    H_GPS = [H_GPS ones(size(H_GPS,1),1)];
end

UB = 0;

end