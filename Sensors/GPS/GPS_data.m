close all
clear all
clc

deg2rad = pi/180;
rad2deg = 180/pi;

c = 2.99792458e8;     % [m/s] speed of light
f_L1 = 1575.42e6;     % [Hz] L1 frequency

OmegaE = 7.2921151467e-5; % [rad/s]

param.earth = struct('c',c,'OmegaE',OmegaE);

test = false;
if test
    load GT_data_test.mat
    load SVvis_GT_test.mat
else
    FP = 'FP4';
    eval(['load ../../GT/GT_data_' FP '.mat'])
    eval(['load SVvis_GT_' FP '.mat'])
end

%% Receiver

rng('shuffle');

GPSclockbias         = 3e5; % [m]
GPSclockbiasDot      = -1e2;  % [m/s]
GPSclockbiasDotDrift = 10;  % [m/s]

GPSphasePSD  = 1.1e-19*c^2; % [m^2/s^2/Hz] - GPS clock phase PSD
GPSfreqPSD   = 4.3e-20*c^2; % [m^2/s^4/Hz] - GPS clock frequency PSD
GPSnoise_rho = 10.5;%10/3;     % [m^2/s^4/Hz] - GPS rho white noise
GPSnoise_vel = 2;%1;        % [m^2/s^4/Hz] - GPS rho white noise

M = 20; % Monte Carlo realizations

T_clock = 0.01;
clock_time = 0:T_clock:tspan(end);
N_clock = length(clock_time);
dt = tspan(2)-tspan(1);
ratioTdt = round(T_clock/dt);

T_GPSTC = 1; %0.20; % [s]
GPSTC_time = 0:T_GPSTC:tspan(end);
N_GPS = length(GPSTC_time);
ratioTPRTclk = round(T_GPSTC/T_clock);

dtSV = SVvis_time(2)-SVvis_time(1);
ratioTdtSV = T_GPSTC/dtSV;
ratioTclkdtSV = round(T_clock/dtSV);

lambda_L1 = c/f_L1;

GPSarm = param.SV.a_Rx_b;

% ### Failure analysis ###
GPSFailureType = '';%'ERRONEOUS'; %'LOSS'
SVIDerror = 24;
type = 'step';%step %ramp
param.failure.a = 120;
param.failure.b = 500;

tFi = 5*60; % s
tFe = 7*60; % min
Finterval = [tFi tFe];

PRerror = PRErrorProfile(GPSTC_time,Finterval,type,param);
% -------

param.GPS = struct('T_GPSTC',T_GPSTC,'GPSphasePSD',GPSphasePSD,'GPSfreqPSD',...
    GPSfreqPSD,'GPSnoise_rho',GPSnoise_rho,'GPSnoise_vel',GPSnoise_vel,...
    'GPSclockbias',GPSclockbias,'GPSclockbiasDot',GPSclockbiasDot,...
    'GPSclockbiasDotDrift',GPSclockbiasDotDrift,'T_clock',T_clock,...
    'lambda_L1',lambda_L1,'GPSarm',GPSarm,'GPSFailureType',GPSFailureType,...
    'Finterval',Finterval);

Xclock = zeros(2,N_clock);

SVPos_eT = cell(1,N_GPS);
SVVel_eT = cell(1,N_GPS);
rho_m = cell(1,N_GPS);
vel_m = cell(1,N_GPS);
URA_m = cell(1,N_GPS);
CN0_m = cell(1,N_GPS);

rhoi_m = cell(1,M);
veli_m = cell(1,M);
GPSTCi_data = cell(1,M);
GPSTCibias_GT = cell(1,M);

rho_f = NaN(N_GPS);

X0_clock = [GPSclockbias GPSclockbiasDot]';
rng('shuffle');
for i=1:M
    
    k=1;
    Xclock(1,k) = X0_clock(1);%*(2*rand-1);
    Xclock(2,k) = X0_clock(2);% + GPSclockbiasDotDrift*(2*rand-1);

    mGPSTC=1;
    nSVs = numel(rho_GT{k});
    rho_m{k} = rho_GT{k} + Xclock(1,1)*ones(nSVs,1) + randn(nSVs,1)*GPSnoise_rho;%<-sqrt(R*T)??
    vel_m{k} = -nu_GT{k}*lambda_L1 +  Xclock(2,1)*ones(nSVs,1) + randn(nSVs,1)*GPSnoise_vel;
    
    for k=2:N_clock
        % ----------
        Xclock(:,k) = GPSclock(Xclock(:,k-1),param);
        % ----------
        if mod(k-1,ratioTPRTclk)==0
            
            mGPSTC = mGPSTC + 1;
            idx_SV = (k-1)*ratioTclkdtSV+1;
            nSVs = numel(rho_GT{idx_SV});
            rho_m{mGPSTC} = rho_GT{idx_SV} + Xclock(1,k)*ones(nSVs,1) +...
                randn(nSVs,1)*GPSnoise_rho;%<-sqrt(R*T)??
            vel_m{mGPSTC} = -nu_GT{idx_SV}*lambda_L1 + Xclock(2,k)*ones(nSVs,1) +...
                randn(nSVs,1)*GPSnoise_vel;
            
            if strcmp(GPSFailureType,'LOSS')
                if mGPSTC>=tFi && mGPSTC<=tFe
                    rho_m{mGPSTC} = [];
                end
            end
            
            if strcmp(GPSFailureType,'ERRONEOUS')
                if mGPSTC>=tFi && mGPSTC<=tFe
                    idx_err = find(SVID{idx_SV}==2);
                    rho_m{mGPSTC}(idx_err) = rho_m{mGPSTC}(idx_err) + PRerror(mGPSTC);
                    rho_f(mGPSTC) = rho_m{mGPSTC}(idx_err);
                end
            end
            
        end
    end
    
    rhoi_m{i} = rho_m;
    veli_m{i} = vel_m;
    
    GPSTCi_data{i} = {rho_m,vel_m};
    GPSTCibias_GT{i} = Xclock(:,:);
    i
end

SV_eT = cell(1,N_GPS);
for k=1:N_GPS
    idx_SV = (k-1)*ratioTdtSV+1;
    % SV_eT{k} = SV_e_GT{idx_SV};
    SVPos_eT{k} = SV_e_GT{idx_SV}(1:3,:);
    SVVel_eT{k} = SV_e_GT{idx_SV}(4:6,:);
    
    URA_m{k} = GPSnoise_rho^2*ones(1,numel(rho_m{k}));
    CN0_m{k} = 40*ones(1,numel(rho_m{k}));
end
GPSTC_aux = {SVID,SVPos_eT,SVVel_eT,URA_m,CN0_m};

paramGPS = param.GPS;
if ~test
    % eval(['save GPS_data_' FP '.mat GPS_time SV_eT rhoi_m veli_m X_clock paramGPS'])
    eval(['save GPSTC_data_' FP '.mat GPSTC_time GPSTCi_data GPSTCibias_GT GPSTC_aux GPSFailureType paramGPS -v7.3'])
end

%% GPS Solution

% sample = fix(rand*M);
sample = 1;

rho = rhoi_m{sample};
vel = veli_m{sample};

R_rec = NaN(3,N_GPS);
V_rec = NaN(3,N_GPS);

GPS_LLA = NaN(3,N_GPS);
GPS_NED = NaN(3,N_GPS);
HDOP = NaN(N_GPS,1);
VDOP = NaN(N_GPS,1);
GDOP = NaN(N_GPS,1);

bias = zeros(N_GPS,1);
b_dot = zeros(N_GPS,1);
visibleSVs = zeros(N_GPS,1);

for k=1:N_GPS
    
    visibleSVs(k) = length(rho{k});
    if visibleSVs(k)>=4
        
        R_SV_eT = SVPos_eT{k}(1:3,:);
        % [R_reck,bias(k)] = GPS_LS(R_SV_eT,rho{k});
        weight = ones(size(rho{k}));
        [R_SV_eR,R_rec(:,k),bias(k)] = GPSsolution(R_SV_eT,rho{k},weight,param);
        % R_SV_eR = R_SV_eT;
        GPS_LLA(:,k) = ECEF2LLA(R_rec(:,k));
        
        [V_rec(:,k),b_dot(k)] = GPSvel_LS(R_SV_eR,GPS_LLA(:,k),SVVel_eT{k},vel{k},weight);
        
        D_ne_GPS = DCM(2,-(GPS_LLA(1,k)+pi/2))*DCM(3,GPS_LLA(2,k));
        GPS_NED(:,k) = D_ne_GPS*(R_rec(:,k)-R_rec(:,1));
        
        [LOS,~] = GPSlinear(R_SV_eR,R_rec(:,k));
        H_til = [-LOS*D_ne_GPS' ones(visibleSVs(k),1)];
        
        R = chol(H_til'*H_til);
        A = R\(R'\eye(4));
        diagA = diag(A);
        HDOP(k) = sqrt(diagA(1)+diagA(2));
        VDOP(k) = sqrt(diagA(3));
        GDOP(k) = sqrt(diagA(1)+diagA(2)+diagA(3)+diagA(4));
        
    else
        GPS_LLA(:,k) = GPS_LLA(:,k-1);
        GPS_NED(:,k) = GPS_NED(:,k-1);
    end
    
end

%% Plot

% sample = fix(rand*M);
% sample = 1;

% ### GPS ###
% LLA
figure
subplot(3,1,1)
hold on
stairs(GPSTC_time,GPS_LLA(1,:)*rad2deg)
plot(tspan,X_GT(:,1)*rad2deg)
xlabel('t [s]')
ylabel('\lambda [deg]')
grid on
legend('GPS','GT')

subplot(3,1,2)
hold on
stairs(GPSTC_time,GPS_LLA(2,:)*rad2deg)
plot(tspan,X_GT(:,2)*rad2deg)
xlabel('t [s]')
ylabel('\Lambda [deg]')
grid on

subplot(3,1,3)
hold on
stairs(GPSTC_time,GPS_LLA(3,:))
plot(tspan,X_GT(:,3))
xlabel('t [s]')
ylabel('altitude [m]')
grid on

% NED
figure
subplot(2,1,1)
hold on
stairs(GPS_NED(2,:),GPS_NED(1,:))
plot(Y_GT(:,5),Y_GT(:,4))
xlabel('East [m]')
ylabel('North [m]')
grid on
legend('GPS','GT')
axis equal

subplot(2,1,2)
hold on
stairs(GPSTC_time,-GPS_NED(3,:))
plot(tspan,-Y_GT(:,6))
xlabel('t [s]')
ylabel('UP [m]')
grid on

% VEL
figure
subplot(3,1,1)
hold on
stairs(GPSTC_time,V_rec(1,:))
plot(tspan,X_GT(:,4))
xlabel('t [s]')
ylabel('V_N [m/s]')
grid on
legend('GPS','GT')

subplot(3,1,2)
hold on
stairs(GPSTC_time,V_rec(2,:))
plot(tspan,X_GT(:,5))
xlabel('t [s]')
ylabel('V_E [m/s]')
grid on

subplot(3,1,3)
hold on
stairs(GPSTC_time,V_rec(3,:))
plot(tspan,X_GT(:,6))
xlabel('t [s]')
ylabel('V_D [m/s]')
grid on

% # Visible SVs | HDOP | VDOP
figure
subplot(3,1,1)
hold on
stairs(GPSTC_time,visibleSVs)
xlabel('t [s]')
ylabel('Number of visible SV')
grid on

subplot(3,1,2)
hold on
stairs(GPSTC_time,HDOP)
xlabel('t [s]')
ylabel('HDOP')
grid on

subplot(3,1,3)
hold on
stairs(GPSTC_time,VDOP)
xlabel('t [s]')
ylabel('VDOP')
grid on

%%
% # bias/clock(1) | clock(2)
figure
subplot(2,1,1)
hold on
stairs(clock_time,Xclock(1,:,sample),'LineWidth',2)
stairs(GPSTC_time,bias,'--','LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$b$ ($m$)','FontSize',14,'interpreter','latex')
legend({'ClockModel','GPSsolution'},'FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex')
grid on

subplot(2,1,2)
hold on
stairs(clock_time,Xclock(2,:,sample),'LineWidth',2)
stairs(GPSTC_time,b_dot,'--','LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$\dot{b}$ ($m/s$)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex')
grid on

% pause
% print -depsc2 FP4clockModel.eps
% close