function [Z_GPSLC,newGPSLC,flagGPSLC] = GPSLC_READ(GPSLC_data,mGPSLC,param)

GPS_LLA = GPSLC_data{1}(:,mGPSLC);
Z_GPSLC = {GPS_LLA};
if ~isempty(GPSLC_data{2})
    GPS_Vel = GPSLC_data{2}(:,mGPSLC);
    Z_GPSLC = {GPS_LLA,GPS_Vel};
end
newGPSLC = true;
flagGPSLC = true;

end