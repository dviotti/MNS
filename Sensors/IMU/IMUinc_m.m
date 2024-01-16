function [IMUinc_data,init] = IMUinc_m(Time,IMU_GT,param)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

% rng('shuffle');

T           = param.IMU.T_IMU;
bias0       = param.IMU.IMUbias;
drift0      = param.IMU.IMUdrift;
% biasstd     = param.IMU.IMUbiasIns;
% driftstd    = param.IMU.IMUdriftIns;
Accstd      = param.IMU.Accstd;
RateGyrostd = param.IMU.RateGyrostd;
tau         = param.IMU.IMUtau;
model       = param.IMU.IMUmodel;

bias_init = randn(3,1)*bias0;%ones(3,1)*bias0;%randn(3,1)*bias0;
drift_init = randn(3,1)*drift0;%ones(3,1)*drift0;%randn(3,1)*drift0;
init = [bias_init;drift_init];

bias = bias_init;
drift = drift_init;

% Time = IMU_data(:,1)';
Acc = IMU_GT(:,1:3)';
AngRate = IMU_GT(:,4:6)';

dt = Time(2) - Time(1);

n=T/dt;
% hN1=n*dt;

N=length(Time);                   % size of continuous-time vector with period step[s]
i_max=fix((N-1)/n)+1;             % length of IMU_inc_time 

switch model
    case 1
        a = 1;
        b = 0;
    case 2
        a = 1;
        b = 0.001;   
    case 3
        a = exp(-dt/tau);
        b = sqrt(1-a^2);
    case 4
        a = exp(-dt/tau);
        b = 0;
    otherwise
        a = 0;
        b = 0;
end

% IMUinc_time = zeros(i_max,1);    % initializes IMU_inc_time  
IMUinc_data = zeros(i_max,12);    % initializes IMU_inc_data

Stor_Acc     = zeros(3,1);                   
Stor_AngRate = zeros(3,1);
% Stor_Bias    = zeros(3,1);                   
% Stor_Drift   = zeros(3,1);
IMUinc_data(1,7:9)=bias';
IMUinc_data(1,10:12)=drift';
for i=2:i_max  % first entry in IMU_inc_data is time integral at initial time = null; begin from i=2
    
    % time integral employs backward trapezoidal rule
    for j=(i-2)*n+2:(i-1)*n+1
        Stor_Acc = Stor_Acc + 0.5*(Acc(:,j)+Acc(:,j-1))*dt;                     % integrates nominal signal
        bias = a*bias + randn(3,1)*bias0*b;                                     % bias dynamic - first-order Gauss- Markov model
        Stor_Acc = Stor_Acc + randn(3,1)*Accstd*sqrt(dt) + bias*dt;             % adds bias and brownian motion with (continuous-time stddev)/n to yield correct sampled-time random walk stdev(hN1)= (continuous-time stddev)*sqrt(hN1)
        Stor_AngRate = Stor_AngRate + 0.5*(AngRate(:,j)+AngRate(:,j-1))*dt;     % integrates nominal signal
        drift = a*drift + randn(3,1)*drift0*b;                                  % drift dynamic - first-order Gauss- Markov model
        Stor_AngRate = Stor_AngRate + randn(3,1)*RateGyrostd*sqrt(dt)+drift*dt; % adds drift and brownian motion with (continuous-time stddev)/n to yield correct sampled-time random walk stdev(hN1)= (continuous-time stddev)*sqrt(hN1)
    end
%     bias = bias*exp(-0.01/tau) + randn(3,1)*bias0*sqrt(1-exp(-2*0.01/tau));
%     drift = drift*exp(-0.01/tau) + randn(3,1)*drift0*sqrt(1-exp(-2*0.01/tau));
    % stores time at sampling period hN1 and accumulated integrals;
    % IMUinc_time(i,1)=Time((i-1)*n+1);
    IMUinc_data(i,1:3)=Stor_Acc';
    IMUinc_data(i,4:6)=Stor_AngRate';
    IMUinc_data(i,7:9)=bias';
    IMUinc_data(i,10:12)=drift';
%     IMUinc_data(i,13:15)=Stor_Bias'./T;
%     IMUinc_data(i,16:18)=Stor_Drift'./T;
    
    % zeroes accumulators at each new round of time integration of n samples from IMU_data
    Stor_Acc     = zeros(3,1);
    Stor_AngRate = zeros(3,1);
%     Stor_Bias    = zeros(3,1);                   
%     Stor_Drift   = zeros(3,1);
end

end

