function [euler_dot,VNED_dot] = FlightProfile4(t)
% FlightProfile4  - Typical Flight with Detour (Obstacle)

deg2rad = pi/180;
% rad2deg = 180/pi;

heading = 45;
heading_rad = 45*deg2rad;

powerup = 20;
takeoff = 10;
yawturn = 3;
transition1 = 20;
straight1 = 40;
straight2 = 20;
climb = 20;
straight3 = 240;
descend = 20;
straight4 = 20;
transition2 = 20;

straight = 10;
straightN = 20;

% Rturn = 1e3;
% turn = pi/2*Rturn/80+0.5;

Ralignment = 2e3;
alignment = heading_rad*Ralignment/58.5;
accalign = 58.5^2/Ralignment+0.1454;
% phialign = atan(accalign/9.7864)*rad2deg;

Rturn = 1.5e3;
turn = pi/2*Rturn/58.5;
accturn = 58.5^2/Rturn+0.148;
% phiturn = atan(accturn/9.7864)*rad2deg;

landing = takeoff;

t1 = powerup;
t2 = t1 + takeoff;
t3 = t2 + 10;
t4 = t3 + yawturn; %3
t5 = t4 + 7;
t6 = t5 + transition1;
t7 = t6 + straight1;
t8 = t7 + alignment;
t9 = t8 + straight2;
t10 = t9 + climb;
t11 = t10 + straight3;
t12 = t11 + descend;
t13 = t12 + straight4;

t14 = t13 + turn;
t15 = t14 + straight;
t16 = t15 + turn;
t17 = t16 + straightN;
t18 = t17 + turn;
t19 = t18 + straight;
t20 = t19 + turn;
t21 = t20 + straightN;

t22 = t21 + transition2;
t23 = t22 + landing;

if t>=t1 && t<t2 % Takeoff
  
    tramp  = 0.5; %[s]
    tinter = 4; %[s]
    takeoffAccel = [0 0 -2]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,2,takeoffAccel,t1,t2,tramp,tinter);
    
elseif t>=t2 && t<t3 % hold
  
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t3 && t<t4 % Heading 45 deg CW in Hover
  
    tramp  = 0.5; %[s]
    tinter = 1; %[s]
    yawrate = heading/(yawturn-tramp);
    turnAglesDot = [0 0 yawrate]'; % [deg/s]
    euler_dot = eulerDotProfile(t,2,turnAglesDot,t3,t4,tramp,tinter);
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t4 && t<t5 % hold
  
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));

elseif t>=t5 && t<t6 % Transition

    tramp  = 0.5; %[s]
    tinter1 = 4; %[s]
    tinter2 = 8; %[s]
    transitionAglesDot = [0 -2 0]'; % [deg/s]
    transitionAccel = [3*cos(heading_rad) 3*sin(heading_rad) -2]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,3,transitionAglesDot,t5,t6,tramp,tinter1,tinter2);
    VNED_dot = VNEDdotProfile(t,5,transitionAccel,t5,t6,tramp,tinter1);

elseif t>=t6 && t<t7 % Straight flight
    
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t7 && t<t8 % Waypoint -> Heading CCW towards North
    
    tramp  = 0.5; %[s]
    tinter = 2; %[s]
    phase = pi/2 - heading_rad;
    yawrate = heading/(alignment-tramp);
    turnAglesDot = [-7 1 -yawrate]'; % [deg/s]
    turnAccel = [accalign/sqrt(2) -accalign/sqrt(2) 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,2,turnAglesDot,t7,t8,tramp,tinter);
    VNED_dot = VNEDdotProfile(t,4,turnAccel,t7,t8,tinter,phase);
    
elseif t>=t8 && t<t9 % Straight flight
    
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t9 && t<t10 % Changing level - Climbing and accelerating
  
    tramp  = 0.5; %[s]
    tinter = 5; %[s]
    climbAglesDot = [0 -0.1 0]'; % [deg/s]
    climbAccel = [1 0 -1.63]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,1,climbAglesDot,t9,t10,tramp);
    VNED_dot = VNEDdotProfile(t,5,climbAccel,t9,t10,tramp,tinter);
    
elseif t>=t10 && t<t11 % Straight flight
    
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));

elseif t>=t11 && t<t12 % Changing level - Descending and desaccelerating
  
    tramp  = 0.5; %[s]
    tinter = 5; %[s]
    climbAglesDot = [0 0.1 0]'; % [deg/s]
    climbAccel = [-1 0 1.63]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,1,climbAglesDot,t11,t12,tramp);
    VNED_dot = VNEDdotProfile(t,5,climbAccel,t11,t12,tramp,tinter);
    
elseif t>=t12 && t<t13 % Straight flight
    
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
% ------

elseif t>=t13 && t<t14 % Coordinated Right Turn N->E
    
    tramp  = 0.5; %[s]
    tinter = 2; %[s]
    phase = 0;
    yawrate = 90/(turn-tramp);
    turnAglesDot = [10 1 yawrate]'; % [deg/s]
    turnAccel = [-accturn/sqrt(2) accturn/sqrt(2) 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,2,turnAglesDot,t13,t14,tramp,tinter);
    VNED_dot = VNEDdotProfile(t,3,turnAccel,t13,t14,tinter,phase);
    
elseif t>=t14 && t<t15 % Straight flight toward East
    
%     leveledAglesDot = [0 0 0]'; % [deg/s]
%     leveledAccel = [0 0 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t15 && t<t16 % Coordinated Left Turn E->N
    
    tramp  = 0.5; %[s]
    tinter = 2; %[s]
    phase = 0;
    yawrate = 90/(turn-tramp);
    turnAglesDot = [-10 1 -yawrate]'; % [deg/s]
    turnAccel = [accturn/sqrt(2) -accturn/sqrt(2) 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,2,turnAglesDot,t15,t16,tramp,tinter);
    VNED_dot = VNEDdotProfile(t,4,turnAccel,t15,t16,tinter,phase);
    
elseif t>=t16 && t<t17 % Straight flight toward North
    
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t17 && t<t18 % Coordinated Left Turn N->W
    
    tramp  = 0.5; %[s]
    tinter = 2; %[s]
    phase = 0;
    yawrate = 90/(turn-tramp);
    turnAglesDot = [-10 1 -yawrate]'; % [deg/s]
    turnAccel = [-accturn/sqrt(2) -accturn/sqrt(2) 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,2,turnAglesDot,t17,t18,tramp,tinter);
    VNED_dot = VNEDdotProfile(t,3,turnAccel,t17,t18,tinter,phase);
    
elseif t>=t18 && t<t19 % Straight flight toward West
    
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t19 && t<t20 % Coordinated Right Turn W->N
    
    tramp  = 0.5; %[s]
    tinter = 2; %[s]
    phase = 0; 
    yawrate = 90/(turn-tramp);
    turnAglesDot = [10 1 yawrate]'; % [deg/s]
    turnAccel = [accturn/sqrt(2) accturn/sqrt(2) 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,2,turnAglesDot,t19,t20,tramp,tinter);
    VNED_dot = VNEDdotProfile(t,4,turnAccel,t19,t20,tinter,phase);
    
elseif t>=t20 && t<t21 % Straight flight toward North
  
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));

% ------
    
elseif t>=t21 && t<t22 % Transition
    
    tramp  = 0.5; %[s]
    tinter1 = 4; %[s]
    tinter2 = 8; %[s]
    transitionAglesDot = [0 -2 0]'; % [deg/s]
    transitionAccel = [-3 0 1]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,3,transitionAglesDot,t21,t22,tramp,tinter2,tinter1); % <-- tinter2,tinter1
    VNED_dot = VNEDdotProfile(t,5,transitionAccel,t21,t22,tramp,tinter1);
    
elseif t>=t22 && t<t23 % Landing
    
    tramp  = 0.5; %[s]
    tinter = 4; %[s]
%     landingAglesDot = [0 0 0]'; % [deg/s]
    takeoffAccel = [0 0 2]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,2,takeoffAccel,t22,t23,tramp,tinter);
    
else
    
    euler_dot = zeros(3,1);
    VNED_dot  = zeros(3,1);
    
end

end