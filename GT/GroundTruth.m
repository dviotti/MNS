close all
clear all
clc

deg2rad = pi/180;
rad2deg = 180/pi;
degph2radps=pi/180/3600;

% Paramiters parameters
mu = 3.986004418e14;      % [m3/s2]
OmegaE = 7.2921151467e-5; % [rad/s]
a = 6378137.0;            % [m]
b = 6356752.3;            % [m]
% f = 1/298.25;           % [] Earth flatness
g0 = 9.780327;            % [m/s2]

R0 = a;
f  = (a-b)/a;

param.earth = struct('mu',mu,'OmegaE',OmegaE,'R0',R0,'f',f,'g0',g0);
c = 2.99792458e8;     % [m/s] speed of light
f_L1 = 1575.42e6;     % [Hz] L1 frequency

dt = 1/400;%1/100;%1/4000;%1/40000;
tf = 750;
tspan = 0:dt:tf;
N = length(tspan);

%% State GT

% ### Flight Profile ###
% FlightProfile0  - Stationary or Straight-and-level Flight - 200s
% FlightProfile1A - 3-axis rotation: pitch-yaw-roll         - 150s
% FlightProfile1B - 3-axis rotation: yaw-pitch-roll         - 150s
% FlightProfile2  - Departure Profile                       - 100s
% FlightProfile3  - Typical Flight                          - 500/600s
% FlightProfile4  - Typical Flight with Detour (Obstacle)   - 750s

cd ./FP
FP = @FlightProfile4;
cd ..

param.state = struct('FP',FP);

% ### attitude GT B2NED ###

epsilon = 0*0.001/3600*deg2rad*sign(randn(3,1));

param.attGT = struct('epsilon',epsilon);

roll0 = 0*deg2rad;
pitch0 = 0*deg2rad;
yaw0 = 0*deg2rad;%15*deg2rad;

D0_bl = DCM(1,roll0)*DCM(2,pitch0)*DCM(3,yaw0);

q0 = DCM2quat(D0_bl);

X0_attGT = zeros(7,1);
X0_attGT(1) = roll0;
X0_attGT(2) = pitch0;
X0_attGT(3) = yaw0;
X0_attGT(4:7) = q0;

% ### navGT ###

lat0 = -(23+12/60)*deg2rad;
long0 = -(45+52/60)*deg2rad;
alt0 = 630;%+370;

VN0 = 0;%80;
VE0 = 0;
VD0 = 0;

% Origin in ECEF
Pos0_ECEF = LLA2ECEF(lat0,long0,alt0);

X0_navGT = zeros(6,1);
X0_navGT(1) = lat0;
X0_navGT(2) = long0;
X0_navGT(3) = alt0;
X0_navGT(4) = VN0;
X0_navGT(5) = VE0;
X0_navGT(6) = VD0;

param.navGT = struct('Pos0_ECEF',Pos0_ECEF);

X0_GT = [
    X0_attGT
    X0_navGT
    ];

tic
[X,Y] = ode4xy(@stateGT,tspan,X0_GT,param);
toc

X_attGT = X(:,1:7);
X_navGT = X(:,8:13);

Y_attGT = Y(:,1:9);
Y_navGT = Y(:,10:21);

% attGT
euler_GT = X_attGT(:,1:3);
q_bl_GT  = X_attGT(:,4:7);

omega_bl_b_GT = Y_attGT(:,1:3);
desalignment = Y_attGT(:,4:6);
euler_dot = Y_attGT(:,7:9);

% navGT
lat_GT  = X_navGT(:,1);
long_GT = X_navGT(:,2);
alt_GT  = X_navGT(:,3);
VN_GT   = X_navGT(:,4);
VE_GT   = X_navGT(:,5);
VD_GT   = X_navGT(:,6);

VN_dot   = Y_navGT(:,4);
VE_dot   = Y_navGT(:,5);
VD_dot   = Y_navGT(:,6);

R_ECEF = Y_navGT(:,7:9);
R_NED  = Y_navGT(:,10:12);

%%

figure
hold on
plot(tspan,omega_bl_b_GT(:,1:3)*rad2deg)
xlabel('t [s]')
ylabel('Angular Velocity [deg/s]')
legend('\omega_{b,x}','\omega_{b,y}','\omega_{b,z}')
grid on

figure
hold on
plot(tspan,wrapToPi(euler_GT(:,1:3))*rad2deg)
xlabel('t [s]')
ylabel('Euler Angles [deg]')
legend('\phi','\theta','\psi')
grid on

figure
hold on
plot(tspan,desalignment(:,1:3)*rad2deg*3600)
xlabel('t [s]')
ylabel('Desalignment [arcsec]')
legend('\psi_N','\psi_E','\psi_D')
grid on

figure
hold on
plot(tspan,euler_dot(:,1:3)*rad2deg)
xlabel('t [s]')
ylabel('Euler Angles Dot [deg]')
legend('d\phidt','d\thetadt','d\psidt')
grid on

%%

% Position
figure
subplot(3,1,1)
plot(tspan,lat_GT*rad2deg)
xlabel('t [s]')
ylabel('\lambda [deg]')
grid on

subplot(3,1,2)
plot(tspan,long_GT*rad2deg)
xlabel('t [s]')
ylabel('\Lambda [deg]')
grid on

subplot(3,1,3)
plot(tspan,alt_GT)
xlabel('t [s]')
ylabel('Height [m]')
grid on

% NED
figure
subplot(2,1,1)
plot(R_NED(:,2),R_NED(:,1))
xlabel('East [m]')
ylabel('North [m]')
grid on
axis equal

subplot(2,1,2)
plot(tspan,-R_NED(:,3))
xlabel('t [s]')
ylabel('HAGL [m]')
grid on

% Velocity
figure
subplot(3,1,1)
plot(tspan,VN_GT)
xlabel('t [s]')
ylabel('V_N [m/s]')
grid on

subplot(3,1,2)
plot(tspan,VE_GT)
xlabel('t [s]')
ylabel('V_E [m/s]')
grid on

subplot(3,1,3)
plot(tspan,VD_GT)
xlabel('t [s]')
ylabel('V_D [m/s]')
grid on

figure
hold on
plot(tspan,VN_dot)
plot(tspan,VE_dot)
plot(tspan,VD_dot)
xlabel('t [s]')
ylabel('Accel [m/s^2]')
legend('dV_N/dt','dV_E/dt','dV_D/dt')
grid on

%% AspGT NED

U_AspNED = [X_navGT Y_navGT]';

AspNED = zeros(3,N);
tic
for i=1:N
    AspNED(:,i) = AspGT(U_AspNED(:,i),param); 
end
Asp_l_GT = AspNED';
toc

% Smoothing
% AspNED_x = smooth(AspNED(1,:),0.005,'lowess');
% AspNED_y = smooth(AspNED(2,:),0.005,'lowess');
% AspNED_z = smooth(AspNED(3,:),0.005,'lowess');
% 
% Asp_l_GT = [AspNED_x AspNED_y AspNED_z];

%%

figure
hold on
plot(tspan,Asp_l_GT(:,1))
plot(tspan,Asp_l_GT(:,2))
plot(tspan,Asp_l_GT(:,3))
xlabel('t [s]')
ylabel('Specific force [m/s^2]')
legend('Asp_N','Asp_E','Asp_D')
grid on

%% OmegaGT and AspGT body

U_IMUGT = [X_navGT(:,1) Y_navGT(:,1:2) X_attGT(:,4:7) omega_bl_b_GT Asp_l_GT]';

IMU = zeros(6,N);

tic
for i=1:N
    IMU(:,i) = IMUGT(U_IMUGT(:,i),param);
end
IMU = IMU';
toc

Asp_b_GT = IMU(:,1:3);
RateGyro_bi_b_GT = IMU(:,4:6);

%%

figure
hold on
plot(tspan,Asp_b_GT(:,1))
plot(tspan,Asp_b_GT(:,2))
plot(tspan,Asp_b_GT(:,3))
xlabel('t [s]')
ylabel('Body specific force [m/s^2]')
legend('N_x','N_y','N_z')
grid on

figure
hold on
plot(tspan,RateGyro_bi_b_GT(:,1)*rad2deg)
plot(tspan,RateGyro_bi_b_GT(:,2)*rad2deg)
plot(tspan,RateGyro_bi_b_GT(:,3)*rad2deg)
xlabel('t [s]')
ylabel('Body rates [deg/s]')
legend('p','q','r')
grid on

%%
U_GT = [Asp_b_GT RateGyro_bi_b_GT];
X_GT = [lat_GT long_GT alt_GT VN_GT VE_GT VD_GT q_bl_GT];
Y_GT = [R_ECEF R_NED euler_GT Asp_l_GT omega_bl_b_GT];

% save GT_data.mat tspan U_GT X_GT Y_GT
save('GT_data_FP4.mat','tspan','U_GT','X_GT','Y_GT','-v7.3')

%%

figure
subplot(3,3,[1,2,3])
hold on
plot(tspan,X_GT(:,4),'LineWidth',2)
plot(tspan,X_GT(:,5),'--','LineWidth',2)
plot(tspan,X_GT(:,6),'-.','LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$V_{NED}$ ($m/s$)','FontSize',14,'interpreter','latex')
legend({'$V_N$','$V_E$','$V_D$'},'FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex')
grid on
xlim([0 750])

subplot(3,3,[4,5,6])
hold on
plot(tspan,wrapToPi(Y_GT(:,7))*rad2deg,'LineWidth',2)
plot(tspan,wrapToPi(Y_GT(:,8))*rad2deg,'--','LineWidth',2)
plot(tspan,wrapToPi(Y_GT(:,9))*rad2deg,'-.','LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$Euler$ ($deg$)','FontSize',14,'interpreter','latex')
legend({'$\phi$','$\theta$','$\psi$'},'FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex')
grid on
xlim([0 750])

% Trajectory
% figure
subplot(3,3,7)
plot(Y_GT(:,5)*1e-3,Y_GT(:,4)*1e-3,'LineWidth',2)
xlabel('$East$ ($km$)','FontSize',14,'interpreter','latex')
ylabel('$North$ ($km$)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex')
grid on
axis equal
% xticks(gca,[-40 -30 -20 -10 -5 0 5 10 20 30 40])

subplot(3,3,[8,9])
plot(tspan,X_GT(:,3),'LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$h$ ($m$)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex')
grid on
xlim([0 750])

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])
