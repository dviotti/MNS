function [Z_ALT,newALT,flagALT] = ALT_READ(ALT_data,mALT,param)

Z_ALT = ALT_data(mALT)*0.3048; %ft2m = 0.3048;
newALT = true;
flagALT = true;

end