function [Xdot,Y] = attitudeGT(t,X,U,param)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

epsilon    = param.attGT.epsilon;

roll = X(1);
pitch = X(2);
yaw = X(3);
q = X(4:7);

% euler_dot = eulerDotProfile(t);
euler_dot = U;

K = [
    1         0           -sin(pitch)
    0  cos(roll) sin(roll)*cos(pitch)
    0 -sin(roll) cos(roll)*cos(pitch)
    ];

omega = K*euler_dot;

omega  = omega + epsilon;

Omega_aug = [
    0     -omega'
    omega -skew(omega)
    ];
q_dot = 1/2*Omega_aug*q;

Xdot = [
    euler_dot
    q_dot
    ];

D_bl = DCM(1,roll)*DCM(2,pitch)*DCM(3,yaw);
D_blc = quat2DCM(1,q); % lc == p
D_lcl = D_blc'*D_bl;

psiN = (D_lcl(2,3)-D_lcl(3,2))/2;
psiE = (D_lcl(3,1)-D_lcl(1,3))/2;
psiD = (D_lcl(1,2)-D_lcl(2,1))/2;

misalignment = [
    psiN
    psiE
    psiD
    ];

Y = [
    omega
    misalignment
    euler_dot
    ];

end

