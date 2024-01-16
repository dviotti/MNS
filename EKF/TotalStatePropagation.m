function [Xout,Yout] = TotalStatePropagation(k,Xin,U,param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

T = param.sim.T;

IMUoutput = param.sensors.IMU.IMUoutput;
IMUtau = param.sensors.IMU.IMUtau;

GPSTC_enabled = param.sensors.GPS.GPSTC_enabled;
ALT_enabled = param.sensors.ALT.ALT_enabled;

% Navegation states
XNAVin = Xin(1:14);
if strcmp(IMUoutput,'increments')
    [XNAVout,YNAVout] = S4S(XNAVin,U,param);
    XNAVout = [XNAVout; zeros(3,1)];
elseif strcmp(IMUoutput,'standard')
    [XNAVoutODE,YNAVoutODE] = ode4xyu(@NAVfull,[0 T],XNAVin,U,param);
    XNAVout = XNAVoutODE(end,:)';
    YNAVout = YNAVoutODE(end,:)';
end

% IMU bias and drift states
XIMU = Xin(15:20);
biasIMU = XIMU(1:3);
driftIMU = XIMU(4:6);
biasIMU = exp(-T/IMUtau)*biasIMU;
driftIMU = exp(-T/IMUtau)*driftIMU;
XIMUout = [
    biasIMU
    driftIMU
    ];

Xout = [
    XNAVout
    XIMUout
    ];
% ALT bias state
if ALT_enabled
    XALT = Xin(end);
    biasALT = XALT;
    XALTout = biasALT;
    Xout = [
        Xout
        XALTout
        ];
end

Yout = YNAVout;
end

