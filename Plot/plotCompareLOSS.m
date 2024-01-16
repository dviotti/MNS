close all
clc

deg2rad = pi/180;
rad2deg = 180/pi;
degph2radps=pi/180/3600;
radps2degph = 1/degph2radps;

g0 = 9.780327;            % [m/s2]
mps22mg = 1e3/g0;

TimeFormat = 1;
switch TimeFormat
    case 1
        TS = 1;
        TF = '($s$)';
    case 2
        TS = 60;
        TF = '($min$)';
end
time_stat = time/TS;

[nKF,N] = size(x_til);
M = 200;%param.config.MCreal;

Covi_x = NaN(nKF,N);
Pi = NaN(nKF,N);
for i=1:N
    Covi_x(:,i) = sqrt(diag(Cov_x(:,:,i)));
    Pi(:,i) = sqrt(diag(P(:,:,i)));
end

nanvec = NaN(1,N);
if nKF==17
    E_x = [E_x; nanvec];
    x_til = [x_til; nanvec];
    Covi_x = [Covi_x; nanvec];
    Pi = [Pi; nanvec];
elseif nKF==16
    E_x = [E_x(1:end-1,:); nanvec; nanvec; E_x(end,:)];
    x_til = [x_til(1:end-1,:); nanvec; nanvec; x_til(end,:)];
    Covi_x = [Covi_x(1:end-1,:); nanvec; nanvec; Covi_x(end,:)];
    Pi = [Pi(1:end-1,:); nanvec; nanvec; Pi(end,:)];
end

Ex_plot = E_x;
x_plot = x_til;
Covx_plot = Covi_x;
P_plot = Pi;

%% Plot

q = 1.96; % 95% double-sided chi-squared hypothesis test
% q = 2.58; % 99% double-sided chi-squared hypothesis test 

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
% LOSS2
limy = {[-3 3],[-3 3],[-3 3],[-0.25 0.25],[-0.25 0.25],[-0.25 0.25],[-5 5],[-5 5],[-5 5],...
    [-1 1],[-1 1],[-1 1],[-2 2],[-2 2],[-2 2],[-10 10],[-1 1],[-0.01 0.01]};

% DR, DV, and Psi
figure
for j=1:6
    subplot(2,3,j)
    hold on
    stairs(time_stat,x_plot(j,:)*convU(j),'Color','#77AC30','LineWidth',2)
    stairs(time_stat,q*P_plot(j,:)*convU(j),'Color','#A2142F','LineWidth',2)
    stairs(time_stat,-q*P_plot(j,:)*convU(j),'Color','#A2142F','LineWidth',2)
    stairs(time_stat,q*Covx_plot(j,:)*convU(j),'Color','#0072BD','LineWidth',2)
    stairs(time_stat,-q*Covx_plot(j,:)*convU(j),'Color','#0072BD','LineWidth',2)
    stairs(time_stat,Ex_plot(j,:)*convU(j),'k','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{j},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    xline([300 420],'--k')
    xlim([0 750])
    xticks([0 250 500 750])
    ylim(limy{j})
end

% set(gcf,'position',[680,558,720,640])
set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.5])

% pause
% print -depsc2 FP4_LOSS1.eps
% print -depsc2 FP4_LOSS2.eps

