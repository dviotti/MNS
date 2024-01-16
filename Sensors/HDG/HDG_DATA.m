close all
clear all
clc

deg2rad = pi/180;
rad2deg = 180/pi;

FP = 'FP4';
eval(['load ../../GT/GT_data_' FP '.mat'])

lat_GT  = X_GT(:,1);
long_GT = X_GT(:,2);
alt_GT  = X_GT(:,3);
% V_GT = X_GT(:,4:6)';
q_GT = X_GT(:,7:10)';

%%
dt = tspan(2) - tspan(1);

T_HDG = 5*0.01;
HDG_time = 0:T_HDG:tspan(end);
N_HDG = length(HDG_time);

stdHDG = 0.5*deg2rad;%10; % (deg)

HDG_DATAtype = 'SIM';
param.HDG = struct('T_HDG',T_HDG,'stdHDG',stdHDG,'HDG_DATAtype',HDG_DATAtype);

M = 200;

% MAGi_data = zeros(6,N_MAG,M);
HDGi_data = cell(M,1);
for i=1:M
   
    HDG = NaN(1,N_HDG);
    euler = NaN(3,N_HDG);
    for k=1:N_HDG
        
        idx = fix(HDG_time(k)/dt) + 1;
        WMM_out = WMM(lat_GT(idx),long_GT(idx),alt_GT(idx));
        D = WMM_out(7);
        
        D_bn = quat2DCM(1,q_GT(:,idx));
        euler(:,k) = DCM2euler(D_bn,'ZYX');

        HDG(k) = euler(3,k) - D + stdHDG*randn;
        
    end
    HDGi_data{i} = HDG*rad2deg;
    i
end

paramHDG = param.HDG;
eval(['save HDG_data_' FP '.mat HDG_time HDGi_data paramHDG'])

%%

% sample = fix(rand*M);
sample = 1;

figure
hold on
plot(HDG_time,HDG*rad2deg)
plot(HDG_time,euler(3,:)*rad2deg,'LineWidth',2)
xlabel('time [s]')
ylabel('Heading (deg)')
grid on
