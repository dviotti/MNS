function [Z_GPSTC,newGPSTC,flagGPSTC] = GPSTC_READ(GPSTC_data,mGPSTC,param)

% GPSvel_enabled = param.GPS.GPSvel_enabled;

rho_m = GPSTC_data{1}{mGPSTC};
vel_m = GPSTC_data{2}{mGPSTC};

SVID_m = GPSTC_data{3}{mGPSTC};
SVPos_eT = GPSTC_data{4}{mGPSTC};
SVVel_eT = GPSTC_data{5}{mGPSTC};

URA_m = GPSTC_data{6}{mGPSTC};
CN0_m = GPSTC_data{7}{mGPSTC};

% vel_m = [];

if ~iscolumn(rho_m)
    rho_m = rho_m';
end
if ~iscolumn(URA_m)
    URA_m = URA_m';
end
if ~iscolumn(CN0_m)
    CN0_m = CN0_m';
end
if ~iscolumn(vel_m)
    vel_m = vel_m';
end

if ~isempty(rho_m)
    if ~isempty(vel_m) %GPSvel_enabled
        Z_GPSTC = {SVID_m,SVPos_eT,rho_m,URA_m,CN0_m,SVVel_eT,vel_m};
    else
        Z_GPSTC = {SVID_m,SVPos_eT,rho_m,URA_m,CN0_m};
    end
    newGPSTC = true;
else
    Z_GPSTC = {};
    newGPSTC = false;
end
flagGPSTC = true;
% flagGPSTC = mGPSTC<29*60 || mGPSTC>31*60;

end