function [euler_dot,VNED_dot] = FlightProfile1A(t)

rot = 1;
hold = 30;

t1 = hold;
t2 = t1 + rot;
t3 = t2 + hold - rot;
t4 = t3 + rot;
t5 = t4 + hold - rot;
t6 = t5 + rot;
t7 = t6 + hold - rot + 30;

if t>=t1 && t<t2 % theta
  
    tramp  = 0; %[s]
    angleDot = [0 45 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,1,angleDot,t1,t2,tramp);
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));

elseif t>=t2 && t<t3 % hold
  
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t3 && t<t4 % psi
    
    tramp  = 0; %[s]
    angleDot = [0 0 45]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,1,angleDot,t3,t4,tramp);
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t4 && t<t5 % hold
    
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t5 && t<t6 % phi
    
    tramp  = 0; %[s]
    angleDot = [45 0 0]'; % [m/s/s]
    euler_dot = eulerDotProfile(t,1,angleDot,t5,t6,tramp);
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
elseif t>=t6 && t<t7 % hold
    
    euler_dot = eulerDotProfile(t,0,zeros(3,1));
    VNED_dot = VNEDdotProfile(t,0,zeros(3,1));
    
else
    
    euler_dot = zeros(3,1);
    VNED_dot  = zeros(3,1);
    
end

end