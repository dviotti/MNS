function [x_hat,P] = EKFprediction(X,U,x_hat,P,Q,param)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

% chi2per = param.config.chi2per;
% Discrete Error-state Kalman Filter

% ### Process Model ### -> Navigation equations in NED
[F,Qd] = EKFmatrix(X,U,Q,param);

% ### Prediction ###
x_hat = F*x_hat;
P = F*P*F' + Qd;
chol(P); % checking if P is positive definite
P = (P+P')/2; % Forcing P to be symmetric

end

function [A,Q] = EKFmatrix(X,U,W,param)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

R0     = param.earth.R0;
OmegaE = param.earth.OmegaE;
f      = param.earth.f;
g0     = param.earth.g0;

T = param.sim.T;

% L = param.sim.L;

% tau = param.sensors.IMU.IMUtau;
% biasIns_on = param.sensors.IMU.biasIns_on;
% % model = param.IMU.sensors.IMUmodel;

tau = param.sim.IMUtau;
biasIns_on = param.sim.biasIns_on;
% model = param.IMU.sensors.IMUmodel;

ALT_enabled = param.sensors.ALT.ALT_enabled;
GPSTC_enabled = param.sensors.GPS.GPSTC_enabled;
% GPStype = param.GPS.GPStype;

lat  = X(1);
% long = X(2);
h    = X(3);
VN   = X(4);
VE   = X(5);
% VD   = X(6);
q_bp = X(7:10);

% dAsp_b = U(1:3);
% Asp_b = dAsp_b/T;
Asp_b = U(1:3);

% Oblate Earth
RE = R0*(1+f*(sin(lat)^2));
RN = R0*(1-f*(2-3*sin(lat)^2));
Re = R0*(1-f*sin(lat)^2);
ge = g0*(1+0.0052884*sin(lat)^2);
% g_NED = g0*(1+0.0052884*sin_lat^2)*(1-2*h/Re);

OmegaE_l = [OmegaE*cos(lat) 0 -OmegaE*sin(lat)]';
rho_l = [VE/(RE+h) -VN/(RN+h) -VE*tan(lat)/(RE+h)]';

Omega1_l = OmegaE_l+rho_l;
Omega2_l = 2*OmegaE_l+rho_l;

q_pb = [q_bp(1); -q_bp(2:4)];
D_pb = quat2DCM(1,q_pb); 
D_cb = D_pb; % At the linearization point: x_hat==0

% -------
% dtheta = [
%     x_hat(2)/(RE+h)
%     -x_hat(1)/(RN+h)
%     -x_hat(2)*tan(lat)/(RE+h)
%     ];
% psi = x_hat(7:9);
% D_pc = eye(3) - skew(psi);
% phi = dtheta + psi;
% D_tp = eye(3) - skew(phi);
% 
% D_cb = D_pc'*D_pb; % Approximation: D_pb = D_cb
% -------

Asp_c = D_cb*Asp_b;
% Asp_t = D_tp*D_pb*Asp_b;

% DR column
A11 = -skew(rho_l);
A21 = diag([-1 -1 2])*ge/Re;
A31 = zeros(3);
A41 = zeros(3);
A51 = zeros(3);

% DV column
A12 = eye(3);
A22 = -skew(Omega2_l);
A32 = zeros(3);
A42 = zeros(3);
A52 = zeros(3);

% psi column
A13 = zeros(3);
A23 = skew(Asp_c);%<----change signal
A33 = -skew(Omega1_l);
A43 = zeros(3);
A53 = zeros(3);

% nabla column
A14 = zeros(3);
A24 = D_pb;
A34 = zeros(3);
A44 = (-1/tau)*eye(3);%zeros(3);
A54 = zeros(3);

% epslon column
A15 = zeros(3);
A25 = zeros(3);
A35 = -D_cb;
A45 = zeros(3);
A55 = (-1/tau)*eye(3);%zeros(3);

Ac = [
    A11 A12 A13 A14 A15
    A21 A22 A23 A24 A25
    A31 A32 A33 A34 A35
    A41 A42 A43 A44 A45
    A51 A52 A53 A54 A55
    ];

% ### Noise Matrix ###

% DR column
B11 = eye(3);
B21 = zeros(3);
B31 = zeros(3);
B41 = zeros(3);
B51 = zeros(3);

% DV column
B12 = zeros(3);
B22 = D_pb;
B32 = zeros(3);
B42 = zeros(3);
B52 = zeros(3);

% psi column
B13 = zeros(3);
B23 = zeros(3);
B33 = D_cb;
B43 = zeros(3);
B53 = zeros(3);

% nabla column
B14 = zeros(3);
B24 = zeros(3);
B34 = zeros(3);
B44 = biasIns_on*eye(3);
B54 = zeros(3);

% epslon column
B15 = zeros(3);
B25 = zeros(3);
B35 = zeros(3);
B45 = zeros(3);
B55 = biasIns_on*eye(3);

% Model with bias instability
Bc = [
    B11 B12 B13 B14 B15
    B21 B22 B23 B24 B25
    B31 B32 B33 B34 B35
    B41 B42 B43 B44 B45
    B51 B52 B53 B54 B55
    ];  

% Model with bias instability
% Bc = [
%     B12 B13 B14 B15
%     B22 B23 B24 B25
%     B32 B33 B34 B35
%     B42 B43 B44 B45
%     B52 B53 B54 B55
%     ];  

% if GPS_enabled && strcmp(GPStype,'tightly')
if GPSTC_enabled
    A_GPS = [
        0 1
        0 0
        ];
    Ac = blkdiag(Ac,A_GPS);
    B_GPS = eye(2);
    Bc = blkdiag(Bc,B_GPS);
end

if ALT_enabled
    A_ALT = 0;
    Ac = blkdiag(Ac,A_ALT);
    B_ALT = 1;
    Bc = blkdiag(Bc,B_ALT);
end

% % ### Descrete State Matrix ###
% Nstates = length(Ac);
% norm_fro_prev = 0;
% A = eye(Nstates);
% norm_fro = norm(A,'fro');
% error_norm = abs(norm_fro-norm_fro_prev);
% k = 1;
% while error_norm>1
% 
%     norm_fro_prev = norm_fro;
%     A = A + (Ac*T)^k/factorial(k);
%     norm_fro = norm(A,'fro');
%     error_norm = abs(norm_fro-norm_fro_prev);
% 
%     k=k+1;
%     if k>5
%         disp('oi')
%         break
%     end
% end
% 
% % ### Descrete Covariance Matrix ###
% QQ = Bc*W*Bc';
% Q = QQ*T + (Ac*QQ+QQ*Ac')*T^2/2;

% Linv = eye(size(L))/L;
% Ac = Linv*Ac*L;
% Bc = Linv*Bc;

n = length(Ac);
Qc = Bc*W*Bc';
Theta = [-Ac Qc;zeros(n) Ac'];
Gamma = expm(Theta*T);
A = Gamma((n+1):2*n,(n+1):2*n)';
Q = A*Gamma(1:n,(n+1):2*n);

end