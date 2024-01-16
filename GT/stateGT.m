function [Xdot,Y] = stateGT(t,X,param)

FP = param.state.FP;

% [euler_dot,VNED_dot] = eval([FP '(t)']);
[euler_dot,VNED_dot] = FP(t);

X_att = X(1:7);
U_att = euler_dot;
[Xdot_att,Y_att] = attitudeGT(t,X_att,U_att,param);

X_nav = X(8:13);
U_nav = VNED_dot;
[Xdot_nav,Y_nav] = navGT(t,X_nav,U_nav,param);

Xdot = [
    Xdot_att
    Xdot_nav
    ];

Y = [
    Y_att
    Y_nav
    ];

end