function [euler_dot,VNED_dot] = FlightProfile2(t)

deg2rad = pi/180;

heading = 45;
heading_rad = 45*deg2rad;

powerup = 20;
takeoff = 10;
yawturn = 3;
forward = 20;

t1 = powerup;
t2 = t1 + takeoff;
t3 = t2 + 10;
t4 = t3 + yawturn; %3
t5 = t4 + 7;
t6 = t5 + forward + 2;

if t>=t1 && t<t2 % Takeoff
  
    tramp  = 0.5; %[s]
    tinter = 4; %[s]
    takeoffAccel = [0 0 -2]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,2,takeoffAccel,t1,t2,tramp,tinter);
    
elseif t>=t2 && t<t3 % Straight flight
  
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

elseif t>=t5 && t<t6 % Forward Hover
  
    tramp  = 0.5; %[s]
    transitionAglesDot = [0 -5 0]'; % [deg/s]
    transitionAccel = [2 2 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,1,transitionAglesDot,t5,t5+2,tramp);
    VNED_dot = VNEDdotProfile(t,1,transitionAccel,t5,t6,2);
    
else
    
    euler_dot = zeros(3,1);
    VNED_dot  = zeros(3,1);
    
end

end