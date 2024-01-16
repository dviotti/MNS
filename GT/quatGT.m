function [q_dot,euler] = quatGT(t,q,omega,param)

Omega_aug = [
    0     -omega'
    omega -skew(omega)
    ];
q_dot = 1/2*Omega_aug*q;

D = Rot(1,q);
% 3-2-3
phi = atan2(D(2,3),-D(1,3));
theta = acos(D(3,3));
psi = atan2(D(3,2),D(3,1));
euler = [phi theta psi];

end