close all
clear all
clc

deg2rad = pi/180;

FP = 'FP4';
% eval(['load ../../GT/GT_data_' FP '.mat'])
eval(['load ../../GT/GT_data_' FP '_0p025.mat'])

T_PVA = 0.2; % 10 ms
PVA_time = 0:T_PVA:tspan(end);
N_PVA = length(PVA_time);

dt = tspan(2) - tspan(1);
ratio = T_PVA/dt;

stdPVApos = 1; % m
stdPVAvel = 0.5; % m/s
stdPVAatt = 1e-3*deg2rad; % rad

param.PVA = struct('T_PVA',T_PVA,'stdPVApos',stdPVApos,...
    'stdPVAvel',stdPVAvel,'stdPVAatt',stdPVAatt);

% tend = 5*60; % s
% kend = fix(tend/T_PVA)+1;
M = 50;

Pos = NaN(3,N_PVA);
Vel = NaN(3,N_PVA);
Att = NaN(3,N_PVA);
PVAi_data = cell(1,M);
for i=1:M
    
   for k=1:N_PVA
    
       idx_GT = (k-1)*ratio+1;
       lat  = X_GT(idx_GT,1);
       long = X_GT(idx_GT,2);
       alt  = X_GT(idx_GT,3);
       D_le = DCM(2,-(lat+pi/2))*DCM(3,long);     
       Pos(:,k) = D_le*LLA2ECEF(lat,long,alt) + stdPVApos*randn(3,1);

%        Pos(:,k) = stdPVAvel*randn(3,1);
       Vel(:,k) = X_GT(idx_GT,4:6)' + stdPVAvel*randn(3,1);
%        Vel(k,:) = [0 0 0];
       
   end
   PVAi_data{i} = {Pos,Vel,Att};
   i
end

paramPVA = param.PVA;
eval(['save PVA_data_' FP '.mat PVA_time PVAi_data paramPVA'])
