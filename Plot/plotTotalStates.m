close all

rad2arcsec = 180/pi*3600;
mps22mg = 1e3/g0;
radps2degph = 180/pi*3600;

ft2m = 0.3048;
m2ft = 1/ft2m;

% sample = fix(rand*M);
sample = M;

TimeFormat = 1;
switch TimeFormat
    case 1
        TS = 1;
        TF = '($s$)';
    case 2
        TS = 60;
        TF = '($min$)';
end

% LLA
figure
subplot(3,1,1)
hold on
% plot(time,X_INS(1,:,sample)*rad2deg)
% stairs(GPS_time,GPS_data3(1,:,sample)*rad2deg)
plot(time/TS,X_hat(1,:)*rad2deg,'LineWidth',2)
plot(tspan/TS,X_GT(:,1)*rad2deg,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\lambda$ (deg)','FontSize',14,'interpreter','latex')
grid on
% legend('GT','INS','GPS','EKF')
% legend('GT','GPS','EKF')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
legend({'EKF','GT'},'FontSize',14,'interpreter','latex')

subplot(3,1,2)
hold on
% plot(time,X_INS(2,:,sample)*rad2deg)
% stairs(GPS_time,GPS_data3(2,:,sample)*rad2deg)
plot(time/TS,X_hat(2,:)*rad2deg,'LineWidth',2)
plot(tspan/TS,X_GT(:,2)*rad2deg,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\Lambda$ (deg)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3)
hold on
% plot(time,X_INS(3,:,sample))
% stairs(GPS_time,GPS_data3(3,:,sample))
plot(time/TS,X_hat(3,:),'LineWidth',2)
plot(tspan/TS,X_GT(:,3),'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$H$ (m)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% Velocity NED
figure
subplot(3,1,1)
hold on
% plot(time,X_INS(4,:,sample))
% stairs(GPS_time,GPS_data3(4,:,sample))
plot(time/TS,X_hat(4,:),'LineWidth',2)
plot(tspan/TS,X_GT(:,4),'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$V_N$ (m/s)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
legend({'EKF','GT'},'FontSize',14,'interpreter','latex')

subplot(3,1,2)
hold on
% plot(time,X_INS(5,:,sample))
% stairs(GPS_time,GPS_data3(5,:,sample))
plot(time/TS,X_hat(5,:),'LineWidth',2)
plot(tspan/TS,X_GT(:,5),'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$V_E$ (m/s)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3)
hold on
% plot(time,X_INS(6,:,sample))
% stairs(GPS_time,GPS_data3(6,:,sample))
plot(time/TS,X_hat(6,:),'LineWidth',2)
plot(tspan/TS,X_GT(:,6),'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$V_D$ (m/s)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% Position NED 
figure
subplot(3,1,[1,2])
hold on
plot(Y_hat(5,:)*1e-3,Y_hat(4,:)*1e-3,'LineWidth',2)
plot(Y_GT(:,5)*1e-3,Y_GT(:,4)*1e-3,'k--','LineWidth',2)
xlabel('$East$ (km)','FontSize',14,'interpreter','latex')
ylabel('$North$ (km)','FontSize',14,'interpreter','latex')
grid on
axis equal
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
legend({'EKF','GT'},'FontSize',14,'interpreter','latex')

subplot(3,1,3)
hold on
plot(time/TS,Y_hat(6,:),'LineWidth',2)
plot(tspan/TS,Y_GT(:,6),'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$Down$ (m)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% Euler angles
figure
subplot(3,1,1)
hold on
plot(time/TS,Y_hat(7,:)*rad2deg,'LineWidth',2)
plot(tspan/TS,Y_GT(:,7)*rad2deg,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\phi$ (deg)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');
legend({'EKF','GT'},'FontSize',14,'interpreter','latex')

subplot(3,1,2)
hold on
plot(time/TS,Y_hat(8,:)*rad2deg,'LineWidth',2)
plot(tspan/TS,Y_GT(:,8)*rad2deg,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\theta$ (deg)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3)
hold on
plot(time/TS,Y_hat(9,:)*rad2deg,'LineWidth',2)
plot(tspan/TS,Y_GT(:,9)*rad2deg,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\psi$ (deg)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

%% IMU 
   
% IMU input
figure
subplot(2,1,1)
plot(IMU_time/TS,IMU_data(1:3,:),'LineWidth',2)
if strcmp(param.sensors.IMU.IMUoutput,'standard')
    ylabel('$f_b (m/s^2)$','FontSize',14,'interpreter','latex')
elseif strcmp(param.sensors.IMU.IMUoutput,'increments')
    ylabel('$\Delta V$ (m/s)','FontSize',14,'interpreter','latex')
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
grid on
legend('X','Y','Z')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(2,1,2)
plot(IMU_time/TS,IMU_data(4:6,:)*rad2deg,'LineWidth',2)
if strcmp(param.sensors.IMU.IMUoutput,'standard')
   ylabel('$\omega$ (deg/s)','FontSize',14,'interpreter','latex')
elseif strcmp(param.sensors.IMU.IMUoutput,'increments')
   ylabel('$\Delta\theta$ (deg)','FontSize',14,'interpreter','latex')
end
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
grid on
legend('X','Y','Z')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% IMU bias
figure
subplot(3,1,1)
hold on
plot(time/TS,X_hat(15,:)*mps22mg,'LineWidth',2)
plot(IMU_time/TS,IMUbias_GT(1,:)*mps22mg,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\nabla$ (mg)','FontSize',14,'interpreter','latex')
grid on
legend({'EKF','GT'},'FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,2)
hold on
plot(time/TS,X_hat(16,:)*mps22mg,'LineWidth',2)
plot(IMU_time/TS,IMUbias_GT(2,:)*mps22mg,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\nabla$ (mg)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3)
hold on
plot(time/TS,X_hat(17,:)*mps22mg,'LineWidth',2)
plot(IMU_time/TS,IMUbias_GT(3,:)*mps22mg,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\nabla$ (mg)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

% IMU drift
figure
subplot(3,1,1)
hold on
plot(time/TS,X_hat(18,:)*radps2degph,'LineWidth',2)
plot(IMU_time/TS,IMUdrift_GT(1,:)*radps2degph,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\epsilon$ (deg/h)','FontSize',14,'interpreter','latex')
grid on
legend({'EKF','GT'},'FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,2)
hold on
plot(time/TS,X_hat(19,:)*radps2degph,'LineWidth',2)
plot(IMU_time/TS,IMUdrift_GT(2,:)*radps2degph,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\epsilon$ (deg/h)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3)
hold on
plot(time/TS,X_hat(20,:)*radps2degph,'LineWidth',2)
plot(IMU_time/TS,IMUdrift_GT(3,:)*radps2degph,'k--','LineWidth',2)
xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
ylabel('$\epsilon$ (deg/h)','FontSize',14,'interpreter','latex')
grid on
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex');

%% Other bias states

% GPS
if GPSTC_enabled && ~ALT_enabled
    figure
    subplot(2,1,1)
    hold on
    plot(time/TS,X_hat(21,:),'.','LineWidth',2)
    plot(time/TS,GPSTCbias_GT(1,:),'k--','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$b_{GPS}$ (m)','FontSize',14,'interpreter','latex')
    grid on
    legend({'EKF','GT'},'FontSize',14,'interpreter','latex')
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,1,2)
    hold on
    plot(time/TS,X_hat(22,:),'.','LineWidth',2)
    plot(time/TS,GPSTCbias_GT(2,:),'k--','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$\dot{b}_{GPS}$ (m/s)','FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

% ALT
if ALT_enabled && ~GPSTC_enabled
    figure
    hold on
    plot(time/TS,X_hat(end,:),'.','LineWidth',2)
    plot(time/TS,ALTbias_GT(1,:),'k--','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$b_{ALT}$ (m)','FontSize',14,'interpreter','latex')
    grid on
    legend({'EKF','GT'},'FontSize',14,'interpreter','latex')
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

% GPS and ALT
if GPSTC_enabled && ALT_enabled
    figure
    subplot(3,1,1)
    hold on
    plot(time/TS,X_hat(21,:),'.','LineWidth',2)
    plot(time/TS,GPSTCbias_GT(1,:),'k--','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$b_{GPS}$ (m)','FontSize',14,'interpreter','latex')
    grid on
    legend({'EKF','GT'},'FontSize',14,'interpreter','latex')
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');

    subplot(3,1,2)
    hold on
    plot(time/TS,X_hat(22,:),'.','LineWidth',2)
    plot(time/TS,GPSTCbias_GT(2,:),'k--','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$\dot{b}_{GPS}$ (m/s)','FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');

    subplot(3,1,3)
    hold on
    plot(time/TS,X_hat(end,:),'.','LineWidth',2)
    plot(time/TS,ALTbias_GT(1,:),'k--','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$b_{ALT}$ (m)','FontSize',14,'interpreter','latex')
    grid on
    legend({'EKF','GT'},'FontSize',14,'interpreter','latex')
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end