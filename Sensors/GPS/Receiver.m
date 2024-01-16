function [rho_m,vel_m,clock] = Receiver(rho_GT,nu_GT,clock,param)

T         = param.GPS.T_GPSTC;
phasePSD  = param.GPS.GPSphasePSD;
freqPSD   = param.GPS.GPSfreqPSD;
noise_rho = param.GPS.GPSnoise_rho;
noise_vel  = param.GPS.GPSnoise_vel;
lambda_L1 = param.GPS.lambda_L1;

% c = param.earth.c;

% -----------------
ctr = clock(1);
ctr_dot = clock(2);

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

clock(1) = ctr;
clock(2) = ctr_dot;
% ------------------

rho_m = rho_GT + ctr + randn*noise_rho;%<-sqrt(R*T)??
vel_m = -nu_GT*lambda_L1 + ctr_dot + randn*noise_vel;

end