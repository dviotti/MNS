close all
clear all
clc

deg2rad = pi/180;
rad2deg = 180/pi;
degph2radps=pi/180/3600;

% Parameters parameters
g0 = 9.780327;            % [m/s2]

FP = 'FP4';
eval(['load ../../GT/GT_data_' FP '.mat'])

%% Run

rng('shuffle');

T_IMU = 1/400; % [s] - Sample time
IMUtype = 'HQTG';%'HQTG' 'NAV'

switch IMUtype
    
    case 'JW'
        
        IMUbias  = 1e-3*g0;         % [m/s^2] - Accelerometer bias
        IMUdrift = 100*degph2radps; % [rad/s] - RateGyro drift
        
        Accstd      = 1e-3*g0;         % [m/s^2] - Accelerometer standard deviation
        RateGyrostd = 15*degph2radps;  % [rad/s] - RateGyro standard deviation
        
    case 'LQTG'
        
        IMUbias  = 1e-3*g0;         % [m/s^2] - Accelerometer bias
        IMUdrift = 10*degph2radps;  % [rad/s] - RateGyro drift
        
        Accstd      = 1e-3*g0;         % [m/s^2] - Accelerometer standard deviation
        RateGyrostd = 15*degph2radps;  % [rad/s] - RateGyro standard deviation
        
    case 'HQTG'
        
        IMUbias  = 0.5*1e-3*g0;   % [m/s^2] - Accelerometer bias
        IMUdrift = 1*degph2radps; % [rad/s] - RateGyro drift
        
        Accstd      = 0.5*1e-3*g0;     % [m/s^2] - Accelerometer standard deviation
        RateGyrostd = 1*degph2radps;  % [rad/s] - RateGyro standard deviation
        
    case 'NAV'
        
        IMUbias  = 0.05*1e-3*g0;    % [m/s^2] - Accelerometer bias
        IMUdrift = 0.1*degph2radps; % [rad/s] - RateGyro drift
        
        Accstd      = 0.05*1e-3*g0;     % [m/s^2] - Accelerometer standard deviation => 2.5e-2 mg/sqrt(Hz) 
        RateGyrostd = 0.1*degph2radps;  % [rad/s] - RateGyro standard deviation => 2.1e-4 deg/h/sqrt(Hz)
        
    case 'OBS' %LEE
        
        IMUbias  = 0.1*1e-3*g0;    % [m/s^2] - Accelerometer bias
        IMUdrift = 0.02*degph2radps; % [rad/s] - RateGyro drift
        
        Accstd      = 0*0.005*1e-3*g0;     % [m/s^2] - Accelerometer standard deviation
        RateGyrostd = 0*0.01*degph2radps;  % [rad/s] - RateGyro standard deviation
        
    otherwise
        
        IMUbias  = 0;    % [m/s^2] - Accelerometer bias
        IMUdrift = 0; % [rad/s] - RateGyro drift
        
        Accstd      = 0;     % [m/s^2] - Accelerometer standard deviation
        RateGyrostd = 0;  % [rad/s] - RateGyro standard deviation
        
end

tau = 1500;%150*10;%3600;%100; % [s] - Time constant

IMUmodel = 1; %1-constant bias; 2- Bias Random Walk; 3- Gauss-Markov; 4- GM no noise;
switch IMUmodel
    case 1
        IMUtau = 1e12;% [s] - Time constant
        biasIns_on = 0;
    case 2
        IMUtau = 1e12;% [s] - Time constant
        biasIns_on = 1;
    case 3
        IMUtau = tau;% [s] - Time constant
        biasIns_on = 1;
    case 4
        IMUtau = tau;% [s] - Time constant
        biasIns_on = 0;
end

IMUarm = [0 0 0]';
IMU_DATAtype = 'SIM';
param.IMU = struct('T_IMU',T_IMU,'IMUtype',IMUtype','IMUbias',IMUbias,...
    'IMUmodel',IMUmodel,'IMUdrift',IMUdrift,'IMUtau',IMUtau,'IMUarm',IMUarm,...
    'Accstd',Accstd,'RateGyrostd',RateGyrostd,'biasIns_on',biasIns_on,...
    'IMU_DATAtype',IMU_DATAtype);

M = 20; % Monte Carlo realizations

IMU_time = 0:T_IMU:tspan(end);
N_IMU = length(IMU_time);
% IMU_data = zeros(N_IMU,7,M);
IMUi_data = NaN(N_IMU,6,M);
IMUibias_GT = NaN(N_IMU,6,M);
initi_data = zeros(6,M);

for i=1:M
    [IMUinc_data,init] = IMUinc_m(tspan,U_GT,param);
    IMUinc_data(1,:) = IMUinc_data(2,:);
    % IMU_data(:,:,i) = [IMU_time' IMUinc_data];
    IMUi_data(:,:,i) = IMUinc_data(:,1:6);
    % IMUibias_GT(:,:,i) = IMUinc_data(:,7:12);
    IMUibias_GT(:,:,i) = IMUinc_data(:,7:12);
    initi_data(:,i) = init;
    i
end

%%
paramIMU = param.IMU;
% save IMU_data.mat IMU_time IMUi_data paramIMU
eval(['save IMU_data_' FP '_' IMUtype '_' num2str(IMUmodel) '.mat IMU_time IMUi_data IMUibias_GT paramIMU -v7.3'])
% save('IMU_data_FP4.mat','IMU_time','IMUi_data','paramIMU','-v7.3')

%% Plot
% sample = fix(rand*M);
T_IMU = IMU_time(2)-IMU_time(1);
sample = 1;

% Specific force
figure
subplot(3,1,1)
hold on
stairs(IMU_time,IMUi_data(:,1,sample)/T_IMU,'LineWidth',2)
plot(tspan,U_GT(:,1),'LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$Asp_{b,x}$ ($m/s^2$)','FontSize',14,'interpreter','latex')
legend({'IMU','GT'},'FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
grid on
xlim([0 750])
xticks([0 250 500 750])

subplot(3,1,2)
hold on
stairs(IMU_time,IMUi_data(:,2,sample)/T_IMU,'LineWidth',2)
plot(tspan,U_GT(:,2),'LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$Asp_{b,y}$ ($m/s^2$)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
grid on
xlim([0 750])
xticks([0 250 500 750])

subplot(3,1,3)
hold on
stairs(IMU_time,IMUi_data(:,3,sample)/T_IMU,'LineWidth',2)
plot(tspan,U_GT(:,3),'LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$Asp_{b,z}$ ($m/s^2$)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
grid on
xlim([0 750])
xticks([0 250 500 750])

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.25 0.25 0.35 0.55])

% pause
% print -depsc2 FP4AspIMU.eps

%%
% Rate angles
figure
subplot(3,1,1)
hold on
stairs(IMU_time,IMUi_data(:,4,sample)/T_IMU*rad2deg,'LineWidth',2)
plot(tspan,U_GT(:,4)*rad2deg,'LineStyle','--','LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$\omega_{b,x}^{BI}$ ($deg/s$)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
legend({'IMU','GT'},'FontSize',14,'interpreter','latex')
grid on
xlim([0 750])
xticks([0 250 500 750])

subplot(3,1,2)
hold on
stairs(IMU_time,IMUi_data(:,5,sample)/T_IMU*rad2deg,'LineWidth',2)
plot(tspan,U_GT(:,5)*rad2deg,'LineStyle','--','LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$\omega_{b,y}^{BI}$ ($deg/s$)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
grid on
xlim([0 750])
xticks([0 250 500 750])

subplot(3,1,3)
hold on
stairs(IMU_time,IMUi_data(:,6,sample)/T_IMU*rad2deg,'LineWidth',2)
plot(tspan,U_GT(:,6)*rad2deg,'LineStyle','--','LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$\omega_{b,z}^{BI}$ ($deg/s$)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
grid on
xlim([0 750])
xticks([0 250 500 750])

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.25 0.25 0.35 0.55])

% pause
% print -depsc2 FP4rateIMU.eps

%%
% Bias and Drift
figure
subplot(2,1,1)
hold on
stairs(IMU_time,IMUibias_GT(:,1,sample)*1e3/g0,'LineWidth',2)
stairs(IMU_time,IMUibias_GT(:,2,sample)*1e3/g0,'LineWidth',2)
stairs(IMU_time,IMUibias_GT(:,3,sample)*1e3/g0,'LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$\nabla$ ($mg$)','FontSize',14,'interpreter','latex')
legend({'$\nabla_x$','$\nabla_y$','$\nabla_z$'},...
    'FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
grid on

subplot(2,1,2)
hold on
stairs(IMU_time,IMUibias_GT(:,4,sample)*rad2deg*3600,'LineWidth',2)
stairs(IMU_time,IMUibias_GT(:,5,sample)*rad2deg*3600,'LineWidth',2)
stairs(IMU_time,IMUibias_GT(:,6,sample)*rad2deg*3600,'LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$\varepsilon$ ($deg/h$)','FontSize',14,'interpreter','latex')
legend({'$\varepsilon_x$','$\varepsilon_x$','$\varepsilon_x$'},...
    'FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
grid on

% pause
% print -depsc2 FP4biasIMU.eps

%%
% Bias
figure
subplot(3,1,1)
hold on
stairs(IMU_time,IMUibias_GT(:,1,sample))
plot(IMU_time,initi_data(1,sample)*ones(1,N_IMU))
xlabel('t [s]')
ylabel('bias_x [m/s^2]')
grid on

subplot(3,1,2)
hold on
stairs(IMU_time,IMUibias_GT(:,2,sample))
plot(IMU_time,initi_data(2,sample)*ones(1,N_IMU))
xlabel('t [s]')
ylabel('bias_y [m/s^2]')
grid on

subplot(3,1,3)
hold on
stairs(IMU_time,IMUibias_GT(:,3,sample))
plot(IMU_time,initi_data(3,sample)*ones(1,N_IMU))
xlabel('t [s]')
ylabel('bias_z [m/s^2]')
grid on

% Drift
figure
subplot(3,1,1)
hold on
stairs(IMU_time,IMUibias_GT(:,4,sample)*rad2deg)
plot(IMU_time,initi_data(4,sample)*ones(1,N_IMU)*rad2deg)
xlabel('t [s]')
ylabel('drift_x [deg/s]')
grid on

subplot(3,1,2)
hold on
stairs(IMU_time,IMUibias_GT(:,5,sample)*rad2deg)
plot(IMU_time,initi_data(5,sample)*ones(1,N_IMU)*rad2deg)
xlabel('t [s]')
ylabel('drift_y [deg/s]')
grid on

subplot(3,1,3)
hold on
stairs(IMU_time,IMUibias_GT(:,6,sample)*rad2deg)
plot(IMU_time,initi_data(6,sample)*ones(1,N_IMU)*rad2deg)
xlabel('t [s]')
ylabel('drift_z [deg/s]')
grid on