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
M = param.config.MCreal;
% M = 200;

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

%% NEES and NIS

q = 1.96; % 95% double-sided chi-squared hypothesis test
% q = 2.58; % 99% double-sided chi-squared hypothesis test 

% NEES
ANEES_plot = ANEES;
dof_ANEES = nKF*M;
limInf_ANEES = 1/2*(-q+sqrt(2*dof_ANEES-1))^2;
limSup_ANEES = 1/2*(q+sqrt(2*dof_ANEES-1))^2;

% NEES_TR
ANEES_TR_plot = ANEES_TR;
dof_ANEES_TR = nKF*M;
limInf_ANEES_TR = 1/2*(-q+sqrt(2*dof_ANEES_TR-1))^2/M;
limSup_ANEES_TR = 1/2*(q+sqrt(2*dof_ANEES_TR-1))^2/M;

% NEES_NI
ANEES_NI_plot = ANEES_NI;
dof_ANEES_NI = 15*M;
limInf_ANEES_NI = 1/2*(-q+sqrt(2*dof_ANEES_NI-1))^2/M;
limSup_ANEES_NI = 1/2*(q+sqrt(2*dof_ANEES_NI-1))^2/M;

% NEES_NAV
ANEES_NAV_plot = ANEES_NAV;
dof_ANEES_NAV = 9*M;
limInf_ANEES_NAV = 1/2*(-q+sqrt(2*dof_ANEES_NAV-1))^2/M;
limSup_ANEES_NAV = 1/2*(q+sqrt(2*dof_ANEES_NAV-1))^2/M;

% NEES_IMU
ANEES_IMU_plot = ANEES_IMU;
dof_ANEES_IMU = 6*M;
limInf_ANEES_IMU = 1/2*(-q+sqrt(2*dof_ANEES_IMU-1))^2/M;
limSup_ANEES_IMU = 1/2*(q+sqrt(2*dof_ANEES_IMU-1))^2/M;

% NEES_GPS
ANEES_GPS_plot = ANEES_GPS;
dof_ANEES_GPS = 2*M;
limInf_ANEES_GPS = 1/2*(-q+sqrt(2*dof_ANEES_GPS-1))^2/M;
limSup_ANEES_GPS = 1/2*(q+sqrt(2*dof_ANEES_GPS-1))^2/M;

% NIS
% idx_NIS = ~isnan(kU(:,end));
% ChiSquare_data = ChiSquareSum(idx,:);
% NIS = sum(ChiSquare_data,2)/M;
ANIS_plot = ANIS;

dof_chi_data = DOF;
dof_NIS = sum(dof_chi_data,2);
NIS_bar = dof_NIS/M;
limInf_NIS = 1/2*(-q+sqrt(2*dof_NIS-1)).^2/M;
limSup_NIS = 1/2*(q+sqrt(2*dof_NIS-1)).^2/M;
% ANIS_plot = ANIS_plot./dof_chi_data;
ANIS_plot = ANIS_plot./NIS_bar;
limInf_NIS_plot = limInf_NIS./NIS_bar;
limSup_NIS_plot = limSup_NIS./NIS_bar;

% dof_NIS_min = min(sum(dof_chi_data,2));
% dof_NIS_max = max(sum(dof_chi_data,2));
% limInf_NIS = 1/2*(-q+sqrt(2*dof_NIS_min-1)).^2/M;
% limSup_NIS = 1/2*(q+sqrt(2*dof_NIS_max-1)).^2/M;

%% Plot

labely = {'$\widetilde{\Delta R}_N$ ($m$)','$\widetilde{\Delta R}_E$ ($m$)','$\widetilde{\Delta R}_D$ ($m$)',...
    '$\widetilde{\Delta V}_N$ ($m/s$)','$\widetilde{\Delta V}_E$ ($m/s$)','$\widetilde{\Delta V}_D$ ($m/s$)',...
    '$\widetilde{\Psi}_N$ ($arcmin$)','$\widetilde{\Psi}_E$ ($arcmin$)','$\widetilde{\Psi}_D$ ($arcmin$)',...
    '$\widetilde{\nabla}_{b,x}$ ($mg$)','$\widetilde{\nabla}_{b,y}$ ($mg$)','$\widetilde{\nabla}_{b,z}$ ($mg$)',...
    '$\widetilde{\varepsilon}_{b,x}$ ($deg/h$)','$\widetilde{\varepsilon}_{b,y}$ ($deg/h$)','$\widetilde{\varepsilon}_{b,z}$ ($deg/h$)',...
    '$\widetilde{b}_{GPS}$ ($m$)','$\widetilde{\dot{b}}_{GPS}$ ($m/s$)','$\widetilde{b}_{ALT}$'};
convU = [1 1 1 1 1 1 rad2deg*60 rad2deg*60 rad2deg*60 mps22mg mps22mg mps22mg radps2degph radps2degph radps2degph 1 1 1];
% GPS / GPS_ALT
% limy = {[-10 10],[-10 10],[-10 10],[-1 1],[-1 1],[-0.5 0.5],[-15 15],[-15 15],[-75 75],...
%     [-5 5],[-5 5],[-2 2],[-10 10],[-10 10],[-10 10],[-10 10],[-1 1],[-0.025 0.025]};
% GPS_ALT_MAG
% limy = {[-10 10],[-10 10],[-10 10],[-1 1],[-1 1],[-0.5 0.5],[-10 10],[-10 10],[-10 10],...
%     [-2 2],[-2 2],[-1 1],[-2 2],[-2 2],[-2 2],[-10 10],[-1 1],[-0.02 0.02]};
% GPS_ALT_MAG_CAM
% limy = {[-5 5],[-5 5],[-5 5],[-0.5 0.5],[-0.5 0.5],[-0.5 0.5],[-5 5],[-5 5],[-5 5],...
%     [-2 2],[-2 2],[-1 1],[-4 4],[-4 4],[-4 4],[-10 10],[-1 1],[-0.01 0.01]};
% GPS_ALT_MAG_CAM_LPS
% limy = {[-4 4],[-4 4],[-2.5 2.5],[-0.25 0.25],[-0.25 0.25],[-0.25 0.25],[-5 5],[-5 5],[-5 5],...
%     [-1 1],[-1 1],[-1 1],[-4 4],[-4 4],[-4 4],[-10 10],[-1 1],[-0.01 0.01]};
% ALT_MAG_CAM_LPS
limy = {[-7.5 7.5],[-7.5 7.5],[-2.5 2.5],[-0.4 0.4],[-0.4 0.4],[-0.25 0.25],[-5 5],[-5 5],[-5 5],...
    [-1 1],[-1 1],[-0.5 0.5],[-4 4],[-4 4],[-4 4],[-10 10],[-1 1],[-0.01 0.01]};

% DR, DV, and Psi
figure
for j=1:9
    subplot(3,3,j)
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
    xlim([0 750])
    xticks([0 250 500 750])
    ylim(limy{j})
end

% set(gcf,'position',[680,558,720,640])
set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])

% pause
% print -depsc2 FP4_GPS1.eps
% print -depsc2 FP4_GPS_ALT1.eps
% print -depsc2 FP4_GPS_ALT_MAG1.eps
% print -depsc2 FP4_GPS_ALT_MAG_CAM1.eps
% print -depsc2 FP4_GPS_ALT_MAG_CAM_LPS1.eps
% print -depsc2 FP4_ALT_MAG_CAM_LPS1.eps

%%
% nabla, epsilon, and bias
figure
for j=1:9
    subplot(3,3,j)
    jj = 9 + j;
    hold on
    stairs(time_stat,x_plot(jj,:)*convU(jj),'Color','#77AC30','LineWidth',2)
    stairs(time_stat,-q*P_plot(jj,:)*convU(jj),'Color','#A2142F','LineWidth',2)
    stairs(time_stat,q*P_plot(jj,:)*convU(jj),'Color','#A2142F','LineWidth',2)
    stairs(time_stat,q*Covx_plot(jj,:)*convU(jj),'Color','#0072BD','LineWidth',2)
    stairs(time_stat,-q*Covx_plot(jj,:)*convU(jj),'Color','#0072BD','LineWidth',2)
    stairs(time_stat,Ex_plot(jj,:)*convU(jj),'k','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{jj},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    xlim([0 750])
    xticks([0 250 500 750])
    ylim(limy{jj})
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])

% pause
% print -depsc2 FP4_GPS2.eps
% print -depsc2 FP4_GPS_ALT2.eps
% print -depsc2 FP4_GPS_ALT_MAG2.eps
% print -depsc2 FP4_GPS_ALT_MAG_CAM2.eps
% print -depsc2 FP4_GPS_ALT_MAG_CAM_LPS2.eps
% print -depsc2 FP4_ALT_MAG_CAM_LPS2.eps


%%
figure
% ANEES TR
subplot(3,2,1)
hold on
plot(time_stat,ANEES_TR_plot,'Color','#0072BD','LineWidth',2)
plot(time_stat,nKF*ones(size(time_stat)),'k','LineWidth',2)
plot(time_stat,limInf_ANEES_TR*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
plot(time_stat,limSup_ANEES_TR*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$ANEES \ TR$','FontSize',14,'interpreter','latex')
% ylim([0 25])
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% ANEES NAV
subplot(3,2,2)
hold on
plot(time_stat,ANEES_NAV_plot,'Color','#0072BD','LineWidth',2)
plot(time_stat,9*ones(size(time_stat)),'k','LineWidth',2)
plot(time_stat,limInf_ANEES_NAV*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
plot(time_stat,limSup_ANEES_NAV*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$ANEES \ NAV$','FontSize',14,'interpreter','latex')
% ylim([0 25])
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% ANEES IMU
subplot(3,2,3)
hold on
plot(time_stat,ANEES_IMU_plot,'Color','#0072BD','LineWidth',2)
plot(time_stat,6*ones(size(time_stat)),'k','LineWidth',2)
plot(time_stat,limInf_ANEES_IMU*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
plot(time_stat,limSup_ANEES_IMU*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$ANEES \ IMU$','FontSize',14,'interpreter','latex')
% ylim([0 25])
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% ANEES GPS
subplot(3,2,4)
hold on
plot(time_stat,ANEES_GPS_plot,'Color','#0072BD','LineWidth',2)
plot(time_stat,2*ones(size(time_stat)),'k','LineWidth',2)
plot(time_stat,limInf_ANEES_GPS*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
plot(time_stat,limSup_ANEES_GPS*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$ANEES \ GPS$','FontSize',14,'interpreter','latex')
% ylim([0 25])
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% ANEES_NI
% figure
subplot(3,2,5)
hold on
plot(time_stat,ANEES_NI_plot,'Color','#0072BD','LineWidth',2)
plot(time_stat,15*ones(size(time_stat)),'k','LineWidth',2)
plot(time_stat,limInf_ANEES_NI*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
plot(time_stat,limSup_ANEES_NI*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$ANEES\ NI$','FontSize',14,'interpreter','latex')
% ylim([0 25])
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% ANIS
% figure
subplot(3,2,6)
hold on
plot(time_stat,ANIS_plot,'.','Color','#0072BD','LineWidth',2)
plot(time_stat,ones(size(time_stat)),'k','LineWidth',2)
plot(time_stat,limInf_NIS_plot,'.','Color','#A2142F','LineWidth',2)
plot(time_stat,limSup_NIS_plot,'.','Color','#A2142F','LineWidth',2)
% plot(time_stat,NIS_bar*ones(size(time_stat)),'k','LineWidth',2)
% plot(time_stat,limInf_NIS*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
% plot(time_stat,limSup_NIS*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$ANIS$','FontSize',14,'interpreter','latex')
% ylim([0 2])
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

%% ANEES 
figure
hold on
plot(time_stat,ANEES_plot,'Color','#0072BD','LineWidth',2)
plot(time_stat,nKF*ones(size(time_stat)),'k','LineWidth',2)
plot(time_stat,limInf_ANEES*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
plot(time_stat,limSup_ANEES*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$ANEES$','FontSize',14,'interpreter','latex')
% ylim([0 25])
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

%%

for i=1:18
    m = min(Covx_plot(i,:));
    mTH = q*m*convU(i);
    t = time_stat(Covx_plot(i,:)==m);
    [t mTH]
end

%%

% % DR
% figure
% % sgtitle('\DeltaR_{NED} - sqrt(Cov+) vs sqrt(P+)')
% subplot(3,1,1)
% hold on
% stairs(time_stat,E_x(1,:),'k','LineWidth',2)
% stairs(time_stat,x_plot(1,:),'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(1,1,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(1,1,:))'),'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(1,1,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(1,1,:))'),'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta R_N$ ($m$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% % legend({'$E(\delta x^+)$','$\delta x^+$','$Cov(\delta x^+)^+$','$\sigma^+$'},'FontSize',14,'interpreter','latex')
% 
% subplot(3,1,2)
% hold on
% stairs(time_stat,E_x(2,:),'k','LineWidth',2)
% stairs(time_stat,x_plot(2,:),'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(2,2,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(2,2,:))'),'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(2,2,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(2,2,:))'),'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta R_E$ ($m$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% subplot(3,1,3)
% hold on
% stairs(time_stat,E_x(3,:),'k','LineWidth',2)
% stairs(time_stat,x_plot(3,:),'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(3,3,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(3,3,:))'),'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(3,3,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(3,3,:))'),'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta R_D$ ($m$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% %%
% % DV
% figure
% % sgtitle('\DeltaR_{NED} - sqrt(Cov+) vs sqrt(P+)')
% subplot(3,1,1)
% hold on
% stairs(time_stat,E_x(4,:),'k','LineWidth',2)
% stairs(time_stat,x_plot(4,:),'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(4,4,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(4,4,:))'),'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(4,4,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(4,4,:))'),'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta V_N$ ($m$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% % legend({'E_x^+','\delta x^+','\sqrt(Cov^+)','\sqrt(P^+)'},'FontSize',14,'interpreter','latex')
% 
% subplot(3,1,2)
% hold on
% stairs(time_stat,E_x(5,:),'k','LineWidth',2)
% stairs(time_stat,x_plot(5,:),'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(5,5,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(5,5,:))'),'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(5,5,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(5,5,:))'),'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta V_E$ ($m$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% subplot(3,1,3)
% hold on
% stairs(time_stat,E_x(6,:),'k','LineWidth',2)
% stairs(time_stat,x_plot(6,:),'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(6,6,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(6,6,:))'),'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(6,6,:))'),'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(6,6,:))'),'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta V_D$ ($m$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% %%
% 
% % psi
% figure
% % sgtitle('\DeltaR_{NED} - sqrt(Cov+) vs sqrt(P+)')
% subplot(3,1,1)
% hold on
% stairs(time_stat,E_x(7,:)*rad2deg/60,'k','LineWidth',2)
% stairs(time_stat,x_plot(7,:)*rad2deg/60,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(7,7,:))')*rad2deg/60,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(7,7,:))')*rad2deg/60,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(7,7,:))')*rad2deg/60,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(7,7,:))')*rad2deg/60,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\psi_N$ ($arcmin$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% % legend({'E_x^+','\delta x^+','\sqrt(Cov^+)','\sqrt(P^+)'},'FontSize',14,'interpreter','latex')
% 
% subplot(3,1,2)
% hold on
% stairs(time_stat,E_x(8,:)*rad2deg/60,'k','LineWidth',2)
% stairs(time_stat,x_plot(8,:)*rad2deg/60,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(8,8,:))')*rad2deg/60,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(8,8,:))')*rad2deg/60,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(8,8,:))')*rad2deg/60,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(8,8,:))')*rad2deg/60,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\psi_E$ ($arcmin$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% subplot(3,1,3)
% hold on
% stairs(time_stat,E_x(9,:)*rad2deg/60,'k','LineWidth',2)
% stairs(time_stat,x_plot(9,:)*rad2deg/60,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(9,9,:))')*rad2deg/60,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(9,9,:))')*rad2deg/60,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(9,9,:))')*rad2deg/60,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(9,9,:))')*rad2deg/60,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\psi_D$ ($arcmin$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% %%
% 
% % Nabla
% figure
% % sgtitle('\DeltaR_{NED} - sqrt(Cov+) vs sqrt(P+)')
% subplot(3,1,1)
% hold on
% stairs(time_stat,E_x(10,:)*mps22mg,'k','LineWidth',2)
% stairs(time_stat,x_plot(10,:)*mps22mg,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(10,10,:))')*mps22mg,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(10,10,:))')*mps22mg,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(10,10,:))')*mps22mg,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(10,10,:))')*mps22mg,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta\nabla_{b,x}$ ($mg$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% % legend({'E_x^+','\delta x^+','\sqrt(Cov^+)','\sqrt(P^+)'},'FontSize',14,'interpreter','latex')
% 
% subplot(3,1,2)
% hold on
% stairs(time_stat,E_x(11,:)*mps22mg,'k','LineWidth',2)
% stairs(time_stat,x_plot(11,:)*mps22mg,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(11,11,:))')*mps22mg,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(11,11,:))')*mps22mg,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(11,11,:))')*mps22mg,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(11,11,:))')*mps22mg,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta\nabla_{b,y}$ ($mg$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% subplot(3,1,3)
% hold on
% stairs(time_stat,E_x(12,:)*mps22mg,'k','LineWidth',2)
% stairs(time_stat,x_plot(12,:)*mps22mg,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(12,12,:))')*mps22mg,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(12,12,:))')*mps22mg,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(12,12,:))')*mps22mg,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(12,12,:))')*mps22mg,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta\nabla_{b,z}$ ($mg$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% %%
% 
% radps2degph = 180/pi*3600;
% 
% % Epsilon
% figure
% % sgtitle('\DeltaR_{NED} - sqrt(Cov+) vs sqrt(P+)')
% subplot(3,1,1)
% hold on
% stairs(time_stat,E_x(13,:)*radps2degph,'k','LineWidth',2)
% stairs(time_stat,x_plot(13,:)*radps2degph,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(13,13,:))')*radps2degph,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(13,13,:))')*radps2degph,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(13,13,:))')*radps2degph,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(13,13,:))')*radps2degph,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta\epsilon_{b,x}$ ($deg/h$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% % legend({'E_x^+','\delta x^+','\sqrt(Cov^+)','\sqrt(P^+)'},'FontSize',14,'interpreter','latex')
% 
% subplot(3,1,2)
% hold on
% stairs(time_stat,E_x(14,:)*radps2degph,'k','LineWidth',2)
% stairs(time_stat,x_plot(14,:)*radps2degph,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(14,14,:))')*radps2degph,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(14,14,:))')*radps2degph,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(14,14,:))')*radps2degph,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(14,14,:))')*radps2degph,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta\epsilon_{b,y}$ ($deg/h$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% subplot(3,1,3)
% hold on
% stairs(time_stat,E_x(15,:)*radps2degph,'k','LineWidth',2)
% stairs(time_stat,x_plot(15,:)*radps2degph,'Color','#77AC30','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(Cov_x(15,15,:))')*radps2degph,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,q*sqrt(squeeze(P_plot(15,15,:))')*radps2degph,'Color','#A2142F','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(Cov_x(15,15,:))')*radps2degph,'Color','#0072BD','LineWidth',2)
% stairs(time_stat,-q*sqrt(squeeze(P_plot(15,15,:))')*radps2degph,'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\Delta\epsilon_{b,z}$ ($deg/h$)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% %% GPS and ALT bias
% 
% if param.sensors.GPS.GPSTC_enabled && param.sensors.ALT.ALT_enabled
%     figure
%     % sgtitle('\DeltaR_{NED} - sqrt(Cov+) vs sqrt(P+)')
%     subplot(3,1,1)
%     hold on
%     stairs(time_stat,E_x(16,:),'k','LineWidth',2)
%     stairs(time_stat,x_plot(16,:),'Color','#77AC30','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(Cov_x(16,16,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(P_plot(16,16,:))'),'Color','#A2142F','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(Cov_x(16,16,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(P_plot(16,16,:))'),'Color','#A2142F','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$b_{GPS}$ ($m$)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
%     % legend({'E_x^+','\delta x^+','\sqrt(Cov^+)','\sqrt(P^+)'},'FontSize',14,'interpreter','latex')
% 
%     subplot(3,1,2)
%     hold on
%     stairs(time_stat,E_x(17,:),'k','LineWidth',2)
%     stairs(time_stat,x_plot(17,:),'Color','#77AC30','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(Cov_x(17,17,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(P_plot(17,17,:))'),'Color','#A2142F','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(Cov_x(17,17,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(P_plot(17,17,:))'),'Color','#A2142F','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$\dot{b}_{GPS}$ ($m/s$)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% 
%     subplot(3,1,3)
%     hold on
%     stairs(time_stat,E_x(end,:),'k','LineWidth',2)
%     stairs(time_stat,x_plot(end,:),'Color','#77AC30','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(Cov_x(end,end,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(P_plot(end,end,:))'),'Color','#A2142F','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(Cov_x(end,end,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(P_plot(end,end,:))'),'Color','#A2142F','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$b_{GPS}$ ($m$)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% end
% 
% if param.sensors.GPS.GPSTC_enabled && ~param.sensors.ALT.ALT_enabled
%     figure
%     % sgtitle('\DeltaR_{NED} - sqrt(Cov+) vs sqrt(P+)')
%     subplot(2,1,1)
%     hold on
%     stairs(time_stat,E_x(16,:),'k','LineWidth',2)
%     stairs(time_stat,x_plot(16,:),'Color','#77AC30','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(Cov_x(16,16,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(P_plot(16,16,:))'),'Color','#A2142F','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(Cov_x(16,16,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(P_plot(16,16,:))'),'Color','#A2142F','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$b_{GPS}$ ($m$)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
%     % legend({'E_x^+','\delta x^+','\sqrt(Cov^+)','\sqrt(P^+)'},'FontSize',14,'interpreter','latex')
% 
%     subplot(2,1,2)
%     hold on
%     stairs(time_stat,E_x(17,:),'k','LineWidth',2)
%     stairs(time_stat,x_plot(17,:),'Color','#77AC30','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(Cov_x(17,17,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(P_plot(17,17,:))'),'Color','#A2142F','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(Cov_x(17,17,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(P_plot(17,17,:))'),'Color','#A2142F','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$\dot{b}_{GPS}$ ($m/s$)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% end
% 
% if param.sensors.ALT.ALT_enabled && ~param.sensors.GPS.GPSTC_enabled
%     figure
%     hold on
%     stairs(time_stat,E_x(end,:),'k','LineWidth',2)
%     stairs(time_stat,x_plot(end,:),'Color','#77AC30','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(Cov_x(end,end,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,q*sqrt(squeeze(P_plot(end,end,:))'),'Color','#A2142F','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(Cov_x(end,end,:))'),'Color','#0072BD','LineWidth',2)
%     stairs(time_stat,-q*sqrt(squeeze(P_plot(end,end,:))'),'Color','#A2142F','LineWidth',2)
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel('$b_{GPS}$ ($m$)','FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% end

%%
% % ANEES 
% figure
% hold on
% plot(time_stat,ANEES_plot,'Color','#0072BD','LineWidth',2)
% plot(time_stat,nKF*ones(size(time_stat)),'k','LineWidth',2)
% plot(time_stat,limInf_ANEES*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
% plot(time_stat,limSup_ANEES*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$ANEES$','FontSize',14,'interpreter','latex')
% % ylim([0 25])
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% % ANIS
% figure
% hold on
% plot(time_stat,ANIS_plot,'.','Color','#0072BD','LineWidth',2)
% plot(time_stat,ones(size(time_stat)),'k','LineWidth',2)
% plot(time_stat,limInf_NIS_plot,'.','Color','#A2142F','LineWidth',2)
% plot(time_stat,limSup_NIS_plot,'.','Color','#A2142F','LineWidth',2)
% % plot(time_stat,NIS_bar*ones(size(time_stat)),'k','LineWidth',2)
% % plot(time_stat,limInf_NIS*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
% % plot(time_stat,limSup_NIS*ones(size(time_stat)),'Color','#A2142F','LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$ANIS$','FontSize',14,'interpreter','latex')
% ylim([0 2])
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
