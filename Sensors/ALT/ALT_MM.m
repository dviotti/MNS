function [z,H,R,aux] =  ALT_MM(pressAlt,X,U,param)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

stdALT = param.sensors.ALT.stdALT;
GPSTC_enabled = param.sensors.GPS.GPSTC_enabled;

stdALTgain = param.sensors.ALT.stdALTgain;

alt = X(3);

bias = X(end);

% H = [zeros(1,2) 1 zeros(1,12)];
H = [zeros(1,2) (1+bias) zeros(1,12)];
if GPSTC_enabled
    H = [H zeros(1,2)];
end

% H = [H 1];
H = [H alt];

z = pressAlt - alt - bias*alt;
% z = pressAlt - alt - bias;

R = (stdALTgain*stdALT)^2;

UB = 0;
aux = {'ALT',UB};

end