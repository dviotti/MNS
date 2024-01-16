function [q,euler,D_bv] = TRIAD(X,Y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

D_bv = Y/X;

% Dopt = D_bv;
% for j=1:10
%    Dopt = (3/2)*Dopt - (1/2)*Dopt*(Dopt.')*Dopt;
% end

phi = atan2(D_bv(2,3),D_bv(3,3))*180/pi;
theta = -asin(D_bv(1,3))*180/pi;
psi = atan2(D_bv(1,2),D_bv(1,1))*180/pi;
euler = [phi theta psi]';

r0 = sqrt((trace(D_bv)+1)/4);
r1 = (D_bv(2,3)-D_bv(3,2))/(4*r0);
r2 = (D_bv(3,1)-D_bv(1,3))/(4*r0);
r3 = (D_bv(1,2)-D_bv(2,1))/(4*r0);
q = [r0 r1 r2 r3]';

% phi = atan2(Dopt(2,3),Dopt(3,3))*180/pi;
% theta = -asin(Dopt(1,3))*180/pi;
% psi = atan2(Dopt(1,2),Dopt(1,1))*180/pi;
% euler = [phi theta psi]';
% 
% r0 = sqrt((trace(Dopt)+1)/4);
% r1 = (Dopt(2,3)-Dopt(3,2))/(4*r0);
% r2 = (Dopt(3,1)-Dopt(1,3))/(4*r0);
% r3 = (Dopt(1,2)-Dopt(2,1))/(4*r0);
% q = [r0 r1 r2 r3]';

end

