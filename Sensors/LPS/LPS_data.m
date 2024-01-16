close all
clear all
clc

deg2rad = pi/180;
rad2deg = 180/pi;

FP = 'FP4';
eval(['load ../../GT/GT_data_' FP '.mat'])

lat_GT  = X_GT(:,1);
long_GT = X_GT(:,2);
alt_GT  = X_GT(:,3);
% V_GT = X_GT(:,4:6)';
q_GT = X_GT(:,7:10)';

%%
dt = tspan(2) - tspan(1);

T_LPS = 20*0.01;
LPS_time = 0:T_LPS:tspan(end);
N_LPS = length(LPS_time);

stdLPS = 0.1;
RVPdistTH = 300;

LPS_DATAtype = 'SIM';
param.LPS = struct('T_LPS',T_LPS,'stdLPS',stdLPS,'RVPdistTH',RVPdistTH,...
    'LPS_DATAtype',LPS_DATAtype);

% Vertiport departure LLA
lat_VP1 = lat_GT(1);
long_VP1 = long_GT(1);
alt_VP1 = alt_GT(1);
RVP1_e = LLA2ECEF(lat_VP1,long_VP1,alt_VP1);
D_VP1_ne = DCM(2,-(lat_VP1+pi/2))*DCM(3,long_VP1);

% Vertiport arival LLA
lat_VP2 = lat_GT(end);
long_VP2 = long_GT(end);
alt_VP2 = alt_GT(end);
RVP2_e = LLA2ECEF(lat_VP2,long_VP2,alt_VP2);
D_VP2_ne = DCM(2,-(lat_VP2+pi/2))*DCM(3,long_VP2);

N_antennas = 4;
Ra_n = zeros(3,N_antennas);
Ra_n(:,1) = 20*[cosd(60) sind(60) -5]';
Ra_n(:,2) = 10*[-cosd(45) sind(45) -7]';
Ra_n(:,3) = 15*[-cosd(45) -sind(45) -4]';
Ra_n(:,4) = 5*[cosd(60) -sind(60) -10]';

Rantenna1_e = zeros(3,N_antennas);
Rantenna2_e = zeros(3,N_antennas);
for i=1:N_antennas
    Rantenna1_e(:,i) = D_VP1_ne'*Ra_n(:,i) + RVP1_e;
    Rantenna2_e(:,i) = D_VP2_ne'*Ra_n(:,i) + RVP2_e;
end
Ra_e = [Rantenna1_e Rantenna2_e];
param.LPS.Ra_e = Ra_e;

M = 20;
LPSi_data = cell(M,1);
for i=1:M
    LPS_dist = NaN(N_LPS,2*N_antennas);
    for k=1:N_LPS
        
        idx = fix(LPS_time(k)/dt) + 1;
        Re = LLA2ECEF(lat_GT(idx),long_GT(idx),alt_GT(idx));

        for j=1:N_antennas
            dist1 = norm(Rantenna1_e(:,j)-Re);
            dist2 = norm(Rantenna2_e(:,j)-Re);
            if dist1<RVPdistTH
                LPS_dist(k,j) = dist1 + stdLPS*randn;
            elseif dist2<RVPdistTH
                LPS_dist(k,j+4) = dist2 + stdLPS*randn;
            end
        end
    end
    LPSi_data{i} = LPS_dist;
    i
end

paramLPS = param.LPS;
eval(['save LPS_data_' FP '.mat LPS_time LPSi_data paramLPS'])

%%

% sample = 1;

figure
plot(LPS_time,LPS_dist)


