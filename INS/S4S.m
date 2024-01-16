function [Xout,Yout] = S4S(Xin,U,param,varargin)

R0     = param.earth.R0;
OmegaE = param.earth.OmegaE;
f      = param.earth.f;
g0     = param.earth.g0;

T         = param.sensors.IMU.T_IMU;
% IMUarm    = param.sensors.IMU.IMUarm;

Kvert     = param.INS.K;
Pos0_ECEF = param.INS.Pos0_ECEF;

lat  = Xin(1);
%     long = Xin(2);
h    = Xin(3);
V_v  = Xin(4:6);
q_bv = Xin(7:10);
csi  = Xin(11);

dV_b_m = U(1:3,:);
dTheta_bi_b_m = U(4:6,:);

alpha = cell(1,4);
beta = cell(1,4);
P = cell(1,4);
for j=1:4
    beta{j} = dV_b_m(:,j);% - skew(dTheta_bi_b_m(:,j))*(skew(dTheta_bi_b_m(:,j))*IMUarm);
    alpha{j} = dTheta_bi_b_m(:,j);
    P{j} = skew(alpha{j});
end

% Altimeter
if ~isempty(varargin)
    ha = varargin{1}; %altimeter
else
    ha = param.INS.Pos0_ECEF(3);
end

sin_lat = sin(lat);
cos_lat = cos(lat);
tan_lat = sin_lat/cos_lat;

% Oblate Earth
RE = R0*(1+f*(sin_lat^2));
RN = R0*(1-f*(2-3*sin_lat^2));
Re = R0*(1-f*sin_lat^2);
g_NED = g0*(1+0.0052884*sin_lat^2)*(1-2*h/Re);

% Step 0
VN = V_v(1);
VE = V_v(2);

rho_v = [
    VE/(RE+h)
    -VN/(RN+h)
    -VE/(RE+h)*tan_lat
    ];
OmegaE_v = [
    OmegaE*cos_lat
    0
    -OmegaE*sin_lat
    ];
omega_vi_v = rho_v + OmegaE_v;

omega_vi_v_norm = norm(omega_vi_v);
omega_vi_v_hat = omega_vi_v/omega_vi_v_norm;
q_vnvo = [
    cos((omega_vi_v_norm*4*T)/2)
    omega_vi_v_hat*sin((omega_vi_v_norm*4*T)/2)
    ];

% Step 1
W_prev = zeros(3,1);
for j=1:4 % Sculling correction
    W = beta{j} - skew(alpha{j})*W_prev + W_prev;
    W = beta{j} - skew(alpha{j})*W + W_prev;
    W_prev = W;
end

dVf_b = W;
q_dVf_b = [0; dVf_b];

% Step 2
dphi = alpha{1} + alpha{2} + alpha{3} + alpha{4} ...% coning correction
    + 2/3*(P{1}*alpha{2} + P{3}*alpha{4}) ...
    + 1/2*(P{1} + P{2})*(alpha{3} + alpha{4}) ...
    +1/30*(P{1} - P{2})*(alpha{3} - alpha{4});

dphi_norm = norm(dphi);
dphi_hat = dphi/dphi_norm;
q_bnbo = [
    cos(dphi_norm/2)
    dphi_hat*sin(dphi_norm/2)
    ];

% Step 3
q_vovn = [q_vnvo(1); -q_vnvo(2:4)];
q_bv = prod_quat(prod_quat(q_vovn,q_bv),q_bnbo);

q_vb = [q_bv(1); -q_bv(2:4)];
q_dVf_v = prod_quat(prod_quat(q_bv,q_dVf_b),q_vb);
dVf_v = q_dVf_v(2:4);

% Step 4
g_v = [0 0 g_NED]';

err = h - ha;
vert_stab = (Kvert(2)*err + csi)*[0 0 1]';
% vert_stab = [0 0 0]'; % -----------------------------------------

dVg_v = (-skew(rho_v+2*OmegaE_v)*V_v + g_v + vert_stab)*4*T ;
dV_v = dVg_v + dVf_v;
V_v = V_v + dV_v;

% ### RK4 ###
hk = 4*T;
xk = [Xin(1:3); Xin(11)];
ui = [V_v; err];
f1 = LLA(xk,           ui, param);
f2 = LLA(xk+0.5*hk*f1, ui, param);
f3 = LLA(xk+0.5*hk*f2, ui, param);
f4 = LLA(xk+hk*f3,     ui, param);
Xk = xk + (hk/6)*(f1 + 2*f2 + 2*f3 + f4);

lat = Xk(1);
long = Xk(2);
h = Xk(3);
csi = Xk(4);

% Xout(1) = lat;
% Xout(2) = long;
% Xout(3) = h;
% Xout(4:6) = V_v;
% Xout(7:10) = q_bv;
% Xout(11) = csi;

Xout = [
    lat
    long
    h
    V_v
    q_bv
    csi
    ];

Pos_ECEF = LLA2ECEF(lat,long,h);
D_ve = DCM(2,-(lat+pi/2))*DCM(3,long);
Pos_NED  = D_ve*(Pos_ECEF-Pos0_ECEF);

D_bv = quat2DCM(1,q_bv); % 3-2-1
euler = DCM2euler(D_bv,'ZYX');

% VN = V_v(1);
% VE = V_v(1);
% VD = V_v(1);
% VG = sqrt(VN^2+VE^2);
% FPA = atan2(-VD,VG);
% TKA = atan2(VE,VN);

dAsp_b_CG = 0;
domega_bi_b = 0;
for i=1:4
dAsp_b_CG = dAsp_b_CG + beta{i};
domega_bi_b = domega_bi_b + alpha{i};
end
Asp_b_CG = dAsp_b_CG*T/4;
omega_bi_b = domega_bi_b*T/4;

omega_bv_b = omega_bi_b - D_bv*omega_vi_v;

Yout = [
    Pos_ECEF
    Pos_NED
    euler
    g_NED
    Asp_b_CG
    omega_bv_b
%     VG
%     FPA
%     TKA
    ];

end

function Xdot = LLA(X,U,param)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

R0        = param.earth.R0;
% OmegaE    = param.Earth.OmegaE;
f         = param.earth.f;
% g0        = param.sim.g0;
% Pos0_ECEF = param.sim.Pos0_ECEF;
K         = param.INS.K;

lat = X(1);
% long = X(2);
h = X(3);
% csi = X(4);

VN = U(1);
VE = U(2);
VD = U(3);
err = U(4);

% Oblate Earth
RE = R0*(1+f*(sin(lat))^2);
RN = R0*(1-f*(2-3*(sin(lat))^2));
% Re = R0*(1-f*(sin(lat))^2);
% g_NED = g0*(1+0.0052884*(sin(lat))^2)*(1-2*h/Re);

% Vertical channel stabilization
k1 = K(1);
% k2 = K(2);
k3 = K(3);
% err = h - ha;
csi_dot = k3*err;

% Position
lat_dot = VN/(RN+h);
long_dot = VE/((RE+h)*cos(lat));
h_dot = -VD - k1*err; %- csi;

Xdot = [
    lat_dot
    long_dot
    h_dot
    csi_dot
    ];

end