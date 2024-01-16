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

T_MAG = 5*0.01;
MAG_time = 0:T_MAG:tspan(end);
N_MAG = length(MAG_time);

stdMAG = 100; % (nT)

param.MAG = struct('T_MAG',T_MAG,'stdMAG',stdMAG);

M = 20;

% MAGi_data = zeros(6,N_MAG,M);
MAGi_data = cell(1,M);
Decl_n = NaN(1,N_MAG);
Decl_b_m = NaN(1,N_MAG);

for i=1:M
   
    MAG_aux = NaN(3,N_MAG);
    for k=1:N_MAG
        
        idx = fix(MAG_time(k)/dt) + 1;
        WMM_out = WMM(lat_GT(idx),long_GT(idx),alt_GT(idx));
        B_n = WMM_out(1:3);
        
        D_bn = quat2DCM(1,q_GT(:,idx));
        Bm_b = D_bn*B_n + stdMAG*randn(3,1);
        
        % MAGi_data(:,k,i) = [Bm_b; B_n];
        MAG_aux(:,k) = Bm_b;
        
        Decl_n(k) = WMM_out(7);
        Decl_b_m(k) = atan2(Bm_b(2),Bm_b(1));
    end
    MAGi_data{i} = MAG_aux;
    i
end

paramMAG = param.MAG;
eval(['save MAG_data_' FP '.mat MAG_time MAGi_data paramMAG'])

%%

% sample = fix(rand*M);
sample = 1;

figure
hold on
plot(MAG_time,Decl_n*rad2deg,'LineWidth',2)
plot(MAG_time,Decl_b_m*rad2deg)
plot(tspan,Y_GT(:,9)*rad2deg)
xlabel('time [s]')
ylabel('Declination [deg]')
grid on
legend('n','m','GT')
