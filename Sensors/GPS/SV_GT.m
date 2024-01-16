close all
clear all
clc

path = pwd;
if ismac
    GPS_path = '/Users/daniel/Documents/_dissertation/code/Sensors/GPS/';
elseif ispc
    GPS_path = 'C:\Users\dviotti\Documents\_dissertation\code\Sensors\GPS';
end
if ~strcmp(path,GPS_path)
    eval(['cd ' GPS_path]);
end

deg2rad = pi/180;
rad2deg = 180/pi;

e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';

% Earth's parameters
mu = 3.986004418e14;      % [m3/s2]
OmegaE = 7.2921151467e-5; % [rad/s]
a = 6378137.0;           % [m]
b = 6356752.3;            % [m]

earth.a = a;
earth.b = b;

OmegaE_e = OmegaE*e3; % Earth angular rate vector

% SV parameters - circular orbit
a_SV = 26560000;   % [m] semi-major axis
T_SV = 43082;      % [s] orbital period
i_SV = 55*deg2rad; % [rad] orbit inclination

R_SV_o = a_SV*e1;            % [m] SV position in So
V_SV_I_o = sqrt(mu/a_SV)*e2; % [m/s] SV inertial velocity in So
udot = 2*pi/T_SV;            % [rad/s] orbital angular rate

OrbitData             % GPS constellation

SimCase = 'BC';

%% Run

switch SimCase
    case 'BC' % ### Birdcage and Tests ###
        Dt = 1; % (s)
        SV_time = 0:Dt:2*T_SV;
    case 'FP' % ### Flight Profiles ###
        % Dt = 1/4000; % 0.25 ms
        Dt = 0.001; % 1 ms
        SV_time = -10:Dt:750;
        % SV_time = -10:Dt:3000;
end
N = length(SV_time);
N_SV = length(Sat_OrbitData);

epoch = 0;

R_SV_e = zeros(N,3,N_SV);
V_SV_E_e = zeros(N,3,N_SV);

tic
for j=1:N_SV
    
    u0 = Sat_OrbitData(j,1);
    glan0 = Sat_OrbitData(j,2);
    
    for i=1:N
        
        te = SV_time(i) + epoch;
    
        % So to Se
        D_oe = DCM(3,u0+udot*te)*DCM(1,i_SV)*DCM(3,glan0-OmegaE*te);
        D_eo = D_oe';
        R_SVj_e = D_eo*R_SV_o;
        R_SV_e(i,:,j) = R_SVj_e';
        V_SVj_E_e = D_eo*V_SV_I_o - skew(OmegaE_e)*R_SVj_e;
        V_SV_E_e(i,:,j) = V_SVj_E_e';
    
    end
    j
end
toc

switch SimCase
    case 'BC' % ### Birdcage and Tests ###
        save SV_GT_BC.mat SV_time R_SV_e V_SV_E_e
        birdcase(R_SV_e,V_SV_E_e,earth);
    case 'FP' % ### Flight Profiles ###
        save('SV_GT.mat','SV_time','R_SV_e','V_SV_E_e','-v7.3')
end

% pause
% print -depsc2 Birdcage.eps

%% Birdcage
function birdcase(R_SV_e,V_SV_E_e,earth)

[~,~,N_SV] = size(R_SV_e);

figure
hold on
ellipsoid(0,0,0,earth.a*1e-3,earth.a*1e-3,earth.b*1e-3)
for j=1:N_SV
    plot3(squeeze(R_SV_e(:,1,j))*1e-3,squeeze(R_SV_e(:,2,j))*1e-3,squeeze(R_SV_e(:,3,j))*1e-3);
    plot3(squeeze(R_SV_e(1,1,j))*1e-3,squeeze(R_SV_e(1,2,j))*1e-3,squeeze(R_SV_e(1,3,j))*1e-3,'ko');
    plot3(squeeze(R_SV_e(end,1,j))*1e-3,squeeze(R_SV_e(end,2,j))*1e-3,squeeze(R_SV_e(end,3,j))*1e-3,'rx');

    plot3([R_SV_e(1,1,j)*1e-3 R_SV_e(1,1,j)*1e-3+V_SV_E_e(1,1,j)],...
        [R_SV_e(1,2,j)*1e-3 R_SV_e(1,2,j)*1e-3+V_SV_E_e(1,2,j)],...
        [R_SV_e(1,3,j)*1e-3 R_SV_e(1,3,j)*1e-3+V_SV_E_e(1,3,j)],'k->');
end
xlabel('X_e (km)')
ylabel('Y_e (km)')
zlabel('Z_e (km)')
grid on
axis equal

end
