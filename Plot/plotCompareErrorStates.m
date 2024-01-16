close all
clear all
clc

g0 = 9.780327;            % [m/s2]

rad2deg = 180/pi;
rad2arcsec = 180/pi*3600;
mps22mg = 1e3/g0;
radps2degph = 180/pi*3600;

ft2m = 0.3048;
m2ft = 1/ft2m;

TestNames = {'SIM_1004_1719_FP4_LQTG1_ALT_GPSTC_MAG_LPS_M1',...
    'SIM_1004_1726_FP4_LQTG1_ALT_GPSTC_CAMLMF_MAG_LPS_M1'};

NC = numel(TestNames);

eval(['load ' TestNames{1} '.mat'])
N = numel(time);

TimeFormat = 1;
switch TimeFormat
    case 1
        TS = 1;
        TF = '($s$)';
    case 2
        TS = 60;
        TF = '($min$)';
end
time_plot = time/TS;

T = param.sim.T;
T_GPSTC = 1;%param.sensors.GPS.T_GPSTC;
T_ALT = 0.04;%param.sensors.ALT.T_ALT;

% states
DRN  = NaN(NC,N);
DRE  = NaN(NC,N);
DRD  = NaN(NC,N);
DVN  = NaN(NC,N);
DVE  = NaN(NC,N);
DVD  = NaN(NC,N);
psiN = NaN(NC,N);
psiE = NaN(NC,N);
psiD = NaN(NC,N);
% bias
biasX      = NaN(NC,N);
biasY      = NaN(NC,N);
biasZ      = NaN(NC,N);
driftX     = NaN(NC,N);
driftY     = NaN(NC,N);
driftZ     = NaN(NC,N);
GPSbias    = NaN(NC,N);
GPSbiasDot = NaN(NC,N);
ALTbias    = NaN(NC,N);

clear time X_hat Y_hat x_hat P x y x_til kU sensorSelected FTList param

sample = 1;

for m=1:NC
    
    eval(['load ' TestNames{m} '.mat'])
    
    ErrorType = 1;
    switch ErrorType
        case 1
            x_plot = squeeze(x(:,:,sample));
        case 2
            x_plot = squeeze(x_til(:,:,sample));
    end
    
    DRN(m,:)  = x_plot(1,:);
    DRE(m,:)  = x_plot(2,:);
    DRD(m,:)  = x_plot(3,:);
    DVN(m,:)  = x_plot(4,:);
    DVE(m,:)  = x_plot(5,:);
    DVD(m,:)  = x_plot(6,:);
    psiN(m,:) = x_plot(7,:);
    psiE(m,:) = x_plot(8,:);
    psiD(m,:) = x_plot(9,:);
    % bias
    biasX(m,:)      = x_plot(10,:);
    biasY(m,:)      = x_plot(11,:);
    biasZ(m,:)      = x_plot(12,:);
    driftX(m,:)     = x_plot(13,:);
    driftY(m,:)     = x_plot(14,:);
    driftZ(m,:)     = x_plot(15,:);
    if param.sensors.GPS.GPSTC_enabled
        GPSbias(m,:)    = x_plot(16,:);
        GPSbiasDot(m,:) = x_plot(17,:);
    end
    if param.sensors.ALT.ALT_enabled
        ALTbias(m,:) = x_plot(end,:);
    end
    
    clear time X_hat Y_hat x_hat P x y x_til kU sensorSelected FTList param
    
end

%% ### ERROR STATES ###

% DR, DV, and Psi
figure
subplot(3,3,1)
hold on
for n=1:NC
    plot(time_plot,DRN(n,:),'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Delta R_N$ (m)','FontSize',14,'interpreter','latex')
grid on
legend('NO CAM','WITH CAM')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,2)
hold on
for n=1:NC
    plot(time_plot,DRE(n,:),'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Delta R_E$ (m)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,3)
hold on
for n=1:NC
    plot(time_plot,DRD(n,:),'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Delta R_D$ (m)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

%DV
subplot(3,3,4)
hold on
for n=1:NC
    plot(time_plot,DVN(n,:),'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Delta V_N$ (m/s)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,5)
hold on
for n=1:NC
    plot(time_plot,DVE(n,:),'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Delta V_E$ (m/s)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,6)
hold on
for n=1:NC
    plot(time_plot,DVD(n,:),'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Delta V_D$ (m/s)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% Psi
subplot(3,3,7)
hold on
for n=1:NC
    plot(time_plot,psiN(n,:)*rad2deg,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\psi_N$ ($deg$)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,8)
hold on
for n=1:NC
    plot(time_plot,psiE(n,:)*rad2deg,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\psi_E$ ($deg$)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,9)
hold on
for n=1:NC
    plot(time_plot,psiD(n,:)*rad2deg,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\psi_D$ ($deg$)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% -----

% IMU
% bias
figure
subplot(3,3,1)
hold on
for n=1:NC
    plot(time_plot,biasX(n,:)*mps22mg,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta\nabla_{b,x}$ (mg)','FontSize',14,'interpreter','latex')
grid on
legend('NO CAM','WITH CAM')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,2)
hold on
for n=1:NC
    plot(time_plot,biasY(n,:)*mps22mg,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta\nabla_{b,y}$ (mg)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,3)
hold on
for n=1:NC
    plot(time_plot,biasZ(n,:)*mps22mg,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta\nabla_{b,z}$ (mg)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% drift
subplot(3,3,4)
hold on
for n=1:NC
    plot(time_plot,driftX(n,:)*radps2degph,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta\epsilon_{b,x}$ (deg/h)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,5)
hold on
for n=1:NC
    plot(time_plot,driftY(n,:)*radps2degph,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta\epsilon_{b,y}$ (deg/h)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,6)
hold on
for n=1:NC
    plot(time_plot,driftZ(n,:)*radps2degph,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta\epsilon_{b,z}$ (deg/h)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% GPS
GPSTC_time_plot = RT(time_plot',T,T_GPSTC,'zoh');
subplot(3,3,7)
hold on
for n=1:NC
    GPSbias_plot = RT(GPSbias(n,:)',T,T_GPSTC,'zoh');
    plot(GPSTC_time_plot,GPSbias_plot,'.','LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta b_{GPS}$ (m)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,3,8)
hold on
for n=1:NC
    GPSbiasDot_plot = RT(GPSbiasDot(m,:)',T,T_GPSTC,'zoh');
    plot(GPSTC_time_plot,GPSbiasDot_plot,'.','LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta \dot{b}_{GPS}$ (m/s)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% ALT
ALT_time_plot = RT(time_plot',T,T_ALT,'zoh');
subplot(3,3,9)
hold on
for n=1:NC
    ALTbias_plot = RT(ALTbias(n,:)',T,T_ALT,'zoh');
    plot(ALT_time_plot,ALTbias_plot,'LineWidth',2)
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta b_{ALT}$ (m)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% % GPS
% if GPSTC_enabled
%     figure
%     subplot(2,1,1)
%     GPSbias = RT(x_plot(16,:)',T,param.sensors.GPS.T_GPSTC,'zoh');
%     plot(GPSTC_time_plot,GPSbias,'.','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$\delta b_{GPS}$ (m)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% 
%     subplot(2,1,2)
%     GPSbiasDot = RT(x_plot(17,:)',T,param.sensors.GPS.T_GPSTC,'zoh');
%     plot(GPSTC_time_plot,GPSbiasDot,'.','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$\delta \dot{b}_{GPS}$ (m/s)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% end
% 
% % ALT
% if ALT_enabled
%     figure
%     subplot(2,1,1)
%     plot(time_plot,x_plot(end,:),'.','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$\delta b_{ALT}$ (m)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% end

