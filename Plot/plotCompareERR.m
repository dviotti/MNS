close all
clear all
clc

deg2rad = pi/180;
rad2deg = 180/pi;
degph2radps=pi/180/3600;
radps2degph = 1/degph2radps;

g0 = 9.780327;            % [m/s2]
mps22mg = 1e3/g0;

% TestNames = {'SIM_1026_2255_FP4_HQTG1_GPSTC_ALT_MAG_CAMLMF_LPS_M1_step',...
%     'SIM_1026_2253_FP4_HQTG1_GPSTC_ALT_MAG_CAMLMF_LPS_M1_stepOFF'};

TestNames = {'SIM_1026_2221_FP4_HQTG1_GPSTC_ALT_MAG_CAMLMF_LPS_M1_ramp',...
    'SIM_1026_2224_FP4_HQTG1_GPSTC_ALT_MAG_CAMLMF_LPS_M1_rampOFF'};

NC = numel(TestNames);

TimeFormat = 1;
switch TimeFormat
    case 1
        TS = 1;
        TF = '($s$)';
    case 2
        TS = 60;
        TF = '($min$)';
end

for m=1:NC
    
    eval(['load ' TestNames{m} '.mat'])  
    
    % M = param.config.MCreal;
    [nKF,N] = size(x_til);

    stddev_sta = NaN(nKF,N);
    stddev_hat = NaN(nKF,N);
    for i=1:N
        stddev_hat(:,i) = sqrt(diag(P(:,:,i)));
    end
    stddev_hat_plot{m} = stddev_hat;
    
%     nanvec = NaN(1,N);
%     if nKF==17
%         E_x = [E_x; nanvec];
%         x_til = [x_til; nanvec];
%         stddev_sta = [stddev_sta; nanvec];
%         stddev_hat = [stddev_hat; nanvec];
%     elseif nKF==16
%         E_x = [E_x(1:end-1,:); nanvec; nanvec; E_x(end,:)];
%         x_til = [x_til(1:end-1,:); nanvec; nanvec; x_til(end,:)];
%         stddev_sta = [stddev_sta(1:end-1,:); nanvec; nanvec; stddev_sta(end,:)];
%         stddev_hat = [stddev_hat(1:end-1,:); nanvec; nanvec; stddev_hat(end,:)];
%     end

    x_plot{m} = x_til;
    
    time_plot = time/TS;
    
    S = whos('-file',[TestNames{m} '.mat']);
    for k = 1:length(S)
        eval(['clear ' S(k).name])
    end
end

%% Plot

q = 1.96; % 95% double-sided chi-squared hypothesis test
% q = 2.58; % 99% double-sided chi-squared hypothesis test 

colors = {'#000000','#0072BD'};
colors2 = {'#A2142F','#77AC30'};
style = {'-','--'};
% labely = {'$\Delta R_N$ ($m$)','$\Delta R_E$ ($m$)','$\Delta R_D$ ($m$)',...
%     '$\Delta V_N$ ($m/s^2$)','$\Delta V_E$ ($m/s^2$)','$\Delta V_D$ ($m/s^2$)',...
%     '$\Psi_N$ ($arcmin$)','$\Psi_E$ ($arcmin$)','$\Psi_D$ ($arcmin$)',...
%     '$\nabla_{b,x}$ ($mg$)','$\nabla_{b,y}$ ($mg$)','$\nabla_{b,z}$ ($mg$)',...
%     '$\varepsilon_{b,x}$ ($deg/h$)','$\varepsilon_{b,y}$ ($deg/h$)','$\varepsilon_{b,z}$ ($deg/h$)',...
%     '$b_{GPS}$ ($m$)','$\dot{b}_{GPS}$ ($m/s$)','$b_{ALT}$'};
labely = {'$\widetilde{\Delta R}_N$ ($m$)','$\widetilde{\Delta R}_E$ ($m$)','$\widetilde{\Delta R}_D$ ($m$)',...
    '$\widetilde{\Delta V}_N$ ($m/s$)','$\widetilde{\Delta V}_E$ ($m/s$)','$\widetilde{\Delta V}_D$ ($m/s$)',...
    '$\widetilde{\Psi}_N$ ($arcmin$)','$\widetilde{\Psi}_E$ ($arcmin$)','$\widetilde{\Psi}_D$ ($arcmin$)',...
    '$\widetilde{\nabla}_{b,x}$ ($mg$)','$\widetilde{\nabla}_{b,y}$ ($mg$)','$\widetilde{\nabla}_{b,z}$ ($mg$)',...
    '$\widetilde{\varepsilon}_{b,x}$ ($deg/h$)','$\widetilde{\varepsilon}_{b,y}$ ($deg/h$)','$\widetilde{\varepsilon}_{b,z}$ ($deg/h$)',...
    '$\widetilde{b}_{GPS}$ ($m$)','$\widetilde{\dot{b}}_{GPS}$ ($m/s$)','$\widetilde{b}_{ALT}$'};
convU = [1 1 1 1 1 1 rad2deg*60 rad2deg*60 rad2deg*60 mps22mg mps22mg mps22mg radps2degph radps2degph radps2degph 1 1 1];

% LOSS1
% limy = {[-75 75],[-75 75],[-25 25],[-1 1],[-1 1],[-0.5 0.5],[-10 10],[-10 10],[-10 10],...
%     [-2 2],[-2 2],[-1 1],[-2 2],[-2 2],[-2 2],[-10 10],[-1 1],[-0.02 0.02]};
% ERR
limy = {[-50 50],[-50 50],[-100 100],[-1.2 1.2],[-1.2 1.2],[-1.2 1.2],[-5 5],[-5 5],[-5 5],...
    [-1 1],[-1 1],[-1 1],[-2 2],[-2 2],[-2 2],[-10 10],[-1 1],[-0.01 0.01]};

% DR, DV, and Psi
figure
for j=1:6
    subplot(2,3,j)
    hold on
    for n=1:NC
        stairs(time_plot,x_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
        stairs(time_plot,q*stddev_hat_plot{n}(j,:)*convU(j),style{n},'Color',colors2{n},'LineWidth',2)
        stairs(time_plot,-q*stddev_hat_plot{n}(j,:)*convU(j),style{n},'Color',colors2{n},'LineWidth',2)
%         stairs(time_plot,q*stddev_hat_plot{n}(j,:)*convU(j),'Color','#A2142F','LineWidth',2)
%         stairs(time_plot,-q*stddev_hat_plot{n}(j,:)*convU(j),'Color','#A2142F','LineWidth',2)
    end
    xline(300,'--')
    xline(420,'--')
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{j},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    xlim([0 750])
    xticks([0 250 500 750])
    ylim(limy{j})
end

% set(gcf,'position',[680,558,720,640])
set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.5])

% pause
% print -depsc2 FP4_ERRstep.eps
% print -depsc2 FP4_ERRramp.eps