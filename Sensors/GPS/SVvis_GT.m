close all
clear all
clc

deg2rad = pi/180;
rad2deg = 180/pi;

% Parameters
c = 2.99792458e8;     % [m/s] speed of light
f_L1 = 1575.42e6;     % [Hz] L1 frequency
OmegaE = 7.2921151467e-5; % [rad/s]

param.earth = struct('c',c,'OmegaE',OmegaE);

% SV parameters - circular orbit
% a_SV = 26560000;   % [m] semi-major axis
T_SV = 43082;      % [s] orbital period
% i_SV = 55*deg2rad; % [rad] orbit inclination

% lambda_L1 = c/f_L1;

% ### Flight Profile ###
FP = 'FP4';

test = false;
if test
    load SV_GT_BC.mat
    location = 'NorthHemisphere';
    switch location
        case 'Equator'
            LLA = [0 0 0]*deg2rad;
        case 'NorthHemisphere'
            LLA = [45 0 0]*deg2rad;
        case 'SouthHemisphere'
            LLA = [-45 0 0]*deg2rad;
        case 'NorthPole'
            LLA = [90 0 0]*deg2rad;
        case 'SouthPole'
            LLA = [-90 0 0]*deg2rad;
    end
    tspan = 0:1:2*T_SV;
    N = numel(tspan);

    lat = LLA(1);
    long = LLA(2);
    h = LLA(3);

    % Attitude
    bank = 0*deg2rad;
    D = DCM(1,bank);
    q_RV = DCM2quat(D);

    RV = [LLA zeros(1,3) q_RV'];
    X_GT = repmat(RV,N,1);

    R_rec_eR = LLA2ECEF(lat,long,h);
    D_leR = DCM(2,-(lat+pi/2))*DCM(3,long);
    R_rec_l = D_leR*R_rec_eR;
    R_rec_l(3) = 0;
    Y_GT = repmat(R_rec_l',N,1);

    T_SVvis = 1; %60; % [s]
else
    load SV_GT.mat
    if ismac
        eval(['load ../../GT/GT_data_' FP '.mat'])
    elseif ispc
        eval(['load ..\..\GT\GT_data_' FP '.mat'])
    end
    T_SVvis = 0.005; % 5 ms
end

SV_dt = SV_time(2) - SV_time(1);
idx0 = find(SV_time==0);
[~,~,N_SV] = size(R_SV_e);

elevation = 5;  % (deg) - 5 deg minimum SV elevation to receive datalink (DO-316)
elevation_rad = elevation*deg2rad;

a_Rx_b = [0 0 0]';

param.SV = struct('c',c,'f_L1',f_L1,'elevation_rad',elevation_rad,...
    'SV_dt',SV_dt,'N_SV',N_SV,'idx0',idx0,'OmegaE',OmegaE,'a_Rx_b',a_Rx_b);

%%

SVvis_time = 0:T_SVvis:tspan(end); % time mark
TimeMark = length(SVvis_time);
dt = tspan(2)-tspan(1);
ratioGT = T_SVvis/dt;

chi_plot = NaN(TimeMark,N_SV);
gamma_plot = NaN(TimeMark,N_SV);
isVisible_plot = zeros(TimeMark,N_SV);

SV_e_GT = cell(TimeMark,1);
rho_GT = cell(TimeMark,1);
nu_GT = cell(TimeMark,1);

SVID = cell(1,TimeMark);

HDOP = NaN(1,TimeMark);
VDOP = NaN(1,TimeMark);

for i=1:TimeMark

    idxGT = (i-1)*ratioGT+1;
    [R_SVvis_e,V_SVvis_E_e,range_vis,nu_vis,az_b,el_b,az_l,el_l,isVisible] = ...
        SVvisible(SVvis_time(i),X_GT(idxGT,:)',U_GT(idxGT,:)',R_SV_e,V_SV_E_e,param);
    
    SV_e_GT{i} = [R_SVvis_e; V_SVvis_E_e];
    rho_GT{i} = range_vis';
    nu_GT{i} = nu_vis';
    
    chi_plot(i,:) = az_l;
    gamma_plot(i,:) = el_l;
    isVisible_plot(i,:) = isVisible;

    SVID{i} = find(isVisible);
    
    az_b_vis = az_b(isVisible==1)';

    el_b_vis = el_b(isVisible==1)';
    
    visibleSVs = sum(isVisible);
    H = [-cos(el_b_vis).*sin(az_b_vis) -cos(el_b_vis).*cos(az_b_vis) ...
        -sin(el_b_vis) ones(visibleSVs,1)];
        
    R = chol(H'*H);
    A = R\(R'\eye(4));
    diagA = diag(A);
    HDOP(i) = sqrt(diagA(1)+diagA(2));
    VDOP(i) = sqrt(diagA(3));

end

if test
    save GT_data_test.mat tspan X_GT Y_GT
    save SVvis_GT_test.mat SVvis_time SVID SV_e_GT rho_GT nu_GT
else
    eval(['save SVvis_GT_' FP '.mat SVvis_time SVID SV_e_GT rho_GT nu_GT param'])
end

%%

% Skyplot
figure
polarplot(0,90,'k.','HandleVisibility','off')
ax = gca;
ax.RTick = [0 elevation 30 60 90];
% ax.RTick = [0 30 60 90];
ax.RDir = 'reverse';
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
hold on
for j=1:N_SV
    SVvis = isVisible_plot(:,j);
    if sum(SVvis)==0
        continue
    end
    chi_vis = chi_plot(SVvis==1,j);
    gamma_vis = gamma_plot(SVvis==1,j);
    polarplot(chi_vis,gamma_vis*rad2deg,'.','LineWidth',2,'DisplayName',char("SV "+string(j)))
    polarplot(chi_vis(1),gamma_vis(1)*rad2deg,'ko','LineWidth',2,'HandleVisibility','off')
    polarplot(chi_vis(end),gamma_vis(end)*rad2deg,'rx','LineWidth',2,'HandleVisibility','off')
end
legend('Location','northeast')
theta = linspace(0,2*pi,180);
rho1 = (90-elevation)*ones(size(theta));
rho2 = 90*ones(size(theta));
polarfill(gca,theta, rho1, rho2, 'black', 0.3)

% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex')

% pause
% print -depsc2 FP4skyplot.eps
% close

figure
hold on
for i=1:N_SV
    plot(SVvis_time, i*isVisible_plot(:,i),'.','LineWidth',2)
end

%%

% Visible SVs
isVisibleTime = sum(isVisible_plot,2)';
figure
subplot(2,1,1)
plot(SVvis_time,isVisibleTime,'LineWidth',2)
grid on
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('Visible SVs','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex')
xlim([0 750])
xticks([0 250 500 750])

subplot(2,1,2)
hold on
plot(SVvis_time,HDOP,'LineWidth',2)
plot(SVvis_time,VDOP,'LineWidth',2)
xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
ylabel('$DOP$','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14);
set(gca,'TickLabelInterpreter','latex')
grid on
legend({'$HDOP$','$VDOP$'},'FontSize',14,'interpreter','latex')
xlim([0 750])
xticks([0 250 500 750])

% subplot(3,1,3)
% plot(SVvis_time,VDOP,'LineWidth',2)
% xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
% ylabel('$VDOP$','FontSize',14,'interpreter','latex')
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex')
% grid on

% pause
% print -depsc2 FP4SVmetrics.eps
% close

% figure
% subplot(2,2,1)
% plot(Y_GT(:,2)*1e-3,Y_GT(:,1)*1e-3,'LineWidth',2)
% xlabel('$East$ ($km$)','FontSize',14,'interpreter','latex')
% ylabel('$North$ ($km$)','FontSize',14,'interpreter','latex')
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex')
% grid on
% axis equal
% 
% subplot(2,2,2)
% plot(tspan,X_GT(:,3),'LineWidth',2)
% xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
% ylabel('$Altitude$ ($m$)','FontSize',14,'interpreter','latex')
% set(gca,'FontSize',14);
% set(gca,'TickLabelInterpreter','latex')
% grid on
% % xlim([0 650])

% pause
% print -depsc2 trjectory.eps