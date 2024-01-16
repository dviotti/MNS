close all
clear all
clc

ft2m = 0.3048;
m2ft = 1/ft2m;

FP = 'FP4';
eval(['load ../../GT/GT_data_' FP '.mat'])

T_ALT = 0.04; % 40 ms
ALT_time = 0:T_ALT:tspan(end);
N_ALT = length(ALT_time);

T_SF = 0.01; % 10 ms
SF_time = 0:T_SF:tspan(end);
% N_SF = length(SF_time);

dt = tspan(2) - tspan(1);
ratio = T_ALT/dt;

stdALT = 1; % [m]

ALT_DATAtype = 'SIM';
param.ALT = struct('T_ALT',T_ALT,'stdALT',stdALT,'ALT_DATAtype',ALT_DATAtype);

M = 200;

ALTi_data = cell(1,M);
ALTibias_GT = cell(1,M);

for i=1:M
   ALT_data = NaN(1,N_ALT);
   bias_ALT = NaN(1,N_ALT);
%    SF = NaN(1,N_ALT);
   SF = 0.025;
   for k=1:N_ALT
    
       idx_GT = (k-1)*ratio+1;
       pressALT = X_GT(idx_GT,3) + stdALT*randn;
%        SF(k) = 0.025;
       bias_ALT(k) = SF*X_GT(idx_GT,3);
       ALT_data(k) = (pressALT+bias_ALT(k))*m2ft;
       
   end
   ALTi_data{i} = ALT_data;
   %ALTibias_GT{i} = SF;%bias_ALT;
   ALTibias_GT{i} = SF*ones(size(SF_time));
   i
end

paramALT = param.ALT;
eval(['save ALT_data_' FP '.mat ALT_time ALTi_data ALTibias_GT paramALT'])

%%

% sample = fix(rand*M);
% sample = 1;

figure
plot(tspan,X_GT(:,3)*m2ft)
hold on
plot(ALT_time,ALT_data)
xlabel('time [s]')
ylabel('altitude [ft]')
grid on
legend('GT','ALT')