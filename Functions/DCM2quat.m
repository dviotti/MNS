function q = DCM2quat(D)

r0 = sqrt((trace(D)+1)/4);
r1 = (D(2,3)-D(3,2))/(4*r0);
r2 = (D(3,1)-D(1,3))/(4*r0);
r3 = (D(1,2)-D(2,1))/(4*r0);
q = [r0 r1 r2 r3]';

end