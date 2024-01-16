function Y = IMUGT(U,param)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

OmegaE    = param.earth.OmegaE;

lat        = U(1);
lat_dot    = U(2);
long_dot   = U(3);
q_bl       = U(4:7);
omega_bl_b = U(8:10);
Asp_l      = U(11:13);

% omega
D_bl = quat2DCM(1,q_bl);
omega_li_l = [
    (OmegaE+long_dot)*cos(lat)
    -lat_dot
    -(OmegaE+long_dot)*sin(lat)
    ];
omega_bi_b = omega_bl_b + D_bl*omega_li_l;

% Asp
% q_lb = [q_bl(1); -q_bl(2:4)];
% Asp_l_q = [0; Asp_l];
% Asp_b_q = prod_quat(prod_quat(q_lb,Asp_l_q),q_bl);
% Asp_b = Asp_b_q(2:4);

Asp_b = D_bl*Asp_l;

Y = [
    Asp_b
    omega_bi_b
    ];

end

