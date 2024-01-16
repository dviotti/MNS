function omega_lb_b = omegaGTe(t)

deg2rad = pi/180;

phi_dot = 5*deg2rad;
theta_dot = 0;
psi_dot = 5*deg2rad;

phi0 = 0;
theta0 = 0;
% psi0 = 45*deg2rad;

phi = phi0+phi_dot*t;
theta = theta0+theta_dot*t;
% psi = psi0+psi_dot*t;

K = [
    1 0         -sin(theta)
    0 cos(phi)  sin(phi)*cos(theta)
    0 -sin(phi) cos(phi)*cos(theta)
    ];
angle_dot = [
    phi_dot
    theta_dot
    psi_dot
    ];

omega_lb_b = K*angle_dot;