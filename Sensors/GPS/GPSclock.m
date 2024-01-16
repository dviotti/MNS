function Xout = GPSclock(Xin,param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

T         = param.GPS.T_clock;
phasePSD  = param.GPS.GPSphasePSD;
freqPSD   = param.GPS.GPSfreqPSD;

% c = param.earth.c;

ctr = Xin(1);
ctr_dot = Xin(2);

varPhase = phasePSD*T + freqPSD*T^3/3;
varFreq = freqPSD*T;

ctr = ctr + T*ctr_dot + randn*sqrt(varPhase);
ctr_dot = ctr_dot + randn*sqrt(varFreq);

% jump = 0.001*c;
% % if abs(ctr-jump)<100
% %     ctr = ctr - jump;
% % end
% if ctr<0
%     ctr = ctr + jump;
% end

Xout(1) = ctr;
Xout(2) = ctr_dot;

end

