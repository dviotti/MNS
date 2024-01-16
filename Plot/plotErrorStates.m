close all

g0 = 9.780327;            % [m/s2]

rad2deg = 180/pi;
rad2arcsec = 180/pi*3600;
mps22mg = 1e3/g0;
radps2degph = 180/pi*3600;

ft2m = 0.3048;
m2ft = 1/ft2m;

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

ErrorType = 2;
switch ErrorType
    case 1
        x_plot = x;
        y_plot = y;
    case 2
        x_plot = x_til;
        y_plot = y_til;
end

% if param.sensors.GPS.GPSTC_enabled
%     GPSTC_time = 0:param.sensors.GPS.T_GPSTC:time(end);
%     GPSbias = RT(x_plot(16,:)',param.sim.T,param.sensors.GPS.T_GPSTC,'zoh');
%     GPSbiasDot = RT(x_plot(17,:)',param.sim.T,param.sensors.GPS.T_GPSTC,'zoh');
% end
% 
% if param.sensors.ALT.ALT_enabled
%     ALT_time = 0:param.sensors.ALT.T_ALT:time(end);
%     ALT_SF = RT(x_plot(end,:)',param.sim.T,param.sensors.ALT.T_ALT,'zoh');
% end

%% ### ERROR STATES ###

% figure
% % plot(timeu/TS,xu_til(1:3,:,sample),'.','LineWidth',2)
% plot(time,x(1:3,:,sample),'LineWidth',2)
% xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
% ylabel('$\Delta R$ (m)','FontSize',14,'interpreter','latex')
% grid on
% legend('N','E','D')
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');

% DR, DV, and Psi
figure
subplot(3,1,1)
plot(time_plot,x_plot(1:3,:),'LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Delta R$ (m)','FontSize',14,'interpreter','latex')
grid on
legend('N','E','D')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,2)
plot(time_plot,x_plot(4:6,:),'LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Delta V$ (m/s)','FontSize',14,'interpreter','latex')
grid on
legend('N','E','D')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3)
plot(time_plot,x_plot(7:9,:)*rad2deg,'LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\psi$ ($deg$)','FontSize',14,'interpreter','latex')
grid on
legend('N','E','D')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% IMU
figure
subplot(2,1,1)
plot(time_plot,x_plot(10:12,:)*mps22mg,'LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta\nabla$ (mg)','FontSize',14,'interpreter','latex')
grid on
legend('X','Y','Z')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(2,1,2)
plot(time_plot,x_plot(13:15,:)*radps2degph,'LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\delta\epsilon$ (deg/h)','FontSize',14,'interpreter','latex')
grid on
legend('X','Y','Z')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% GPS
if param.sensors.GPS.GPSTC_enabled && ~param.sensors.ALT.ALT_enabled
    figure
    subplot(2,1,1)
%     plot(GPSTC_time,GPSbias,'.','LineWidth',2)
    plot(time_plot,x_plot(16,:),'LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$\delta b_{GPS}$ (m)','FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,1,2)
%     plot(GPSTC_time,GPSbiasDot,'.','LineWidth',2)
    plot(time_plot,x_plot(17,:),'LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$\delta \dot{b}_{GPS}$ (m/s)','FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

% ALT
if param.sensors.ALT.ALT_enabled && ~param.sensors.GPS.GPSTC_enabled
    figure
%     plot(ALT_time,ALT_SF,'.','LineWidth',2)
    plot(time_plot,x_plot(end,:),'LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$\delta b_{ALT}$ (m)','FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

% GPS
if param.sensors.GPS.GPSTC_enabled && param.sensors.ALT.ALT_enabled
    figure
    subplot(3,1,1)
%     plot(GPSTC_time,GPSbias,'.','LineWidth',2)
    plot(time_plot,x_plot(16,:),'LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$\delta b_{GPS}$ (m)','FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');

    subplot(3,1,2)
%     plot(GPSTC_time,GPSbiasDot,'.','LineWidth',2)
    plot(time_plot,x_plot(17,:),'LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$\delta \dot{b}_{GPS}$ (m/s)','FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');

    subplot(3,1,3)
%     plot(ALT_time,ALT_SF,'.','LineWidth',2)
    plot(time_plot,x_plot(end,:),'LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$\delta b_{ALT}$ (m)','FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

%% Pinson Angles

% figure
% subplot(3,1,1)
% hold on
% plot(time_plot,y_plot(1:3,:)*rad2deg,'LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\phi$ (deg)','FontSize',14,'interpreter','latex')
% grid on
% legend({'X','Y','Z'},'FontSize',14,'interpreter','latex')
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% subplot(3,1,2)
% hold on
% plot(time_plot,y_plot(4:6,:)*rad2deg,'LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\psi$ (deg)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
% 
% subplot(3,1,3)
% hold on
% plot(time_plot,y_plot(7:9,:)*rad2deg,'LineWidth',2)
% xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
% ylabel('$\delta\theta$ (deg)','FontSize',14,'interpreter','latex')
% grid on
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex');
