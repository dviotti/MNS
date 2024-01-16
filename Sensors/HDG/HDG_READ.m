function [Z_HDG,newHDG,flagHDG] = HDG_READ(HDG_data,mHDG,param)

Z_HDG = HDG_data(mHDG)*pi/180;
newHDG = true;
flagHDG = true;

end