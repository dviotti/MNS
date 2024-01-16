function euler_dot = eulerDotProfile(t,profile,varargin)

deg2rad = pi/180;

angularRates = varargin{1};
rollrate  = angularRates(1)*deg2rad;
pitchrate = angularRates(2)*deg2rad;
yawrate   = angularRates(3)*deg2rad;

switch profile
    
    case 1 % Angle Change
        
        tinit = varargin{2};
        tend  = varargin{3};
        tramp = varargin{4};
        
        if t>=tinit && t<(tinit+tramp)
            roll_dot = rollrate*(t-tinit)/tramp;
            pitch_dot = pitchrate*(t-tinit)/tramp;
            yaw_dot = yawrate*(t-tinit)/tramp;
        elseif t>=(tinit+tramp) && t<(tend-tramp)
            roll_dot = rollrate;
            pitch_dot = pitchrate;
            yaw_dot = yawrate;
        elseif t>=(tend-tramp) && t<tend
            roll_dot = rollrate*(tend-t)/tramp;
            pitch_dot = pitchrate*(tend-t)/tramp;
            yaw_dot = yawrate*(tend-t)/tramp;
        else
            roll_dot = 0;
            pitch_dot = 0;
            yaw_dot = 0;
        end
        
    case 2 % Coordinated Turn
        
        tinit = varargin{2};
        tend  = varargin{3};
        tramp = varargin{4};
        tinter = varargin{5};
        
        if t>=tinit && t<(tinit+tramp)
            roll_dot = rollrate*(t-tinit)/tramp;
            pitch_dot = pitchrate*(t-tinit)/tramp;
            yaw_dot = yawrate*(t-tinit)/tramp;
        elseif t>=(tinit+tramp) && t<(tinit+tinter-tramp)
            roll_dot = rollrate;
            pitch_dot = pitchrate;
            yaw_dot = yawrate;
        elseif t>=(tinit+tinter-tramp) && t<(tinit+tinter)
            roll_dot = rollrate*((tinit+tinter)-t)/tramp;
            pitch_dot = pitchrate*((tinit+tinter)-t)/tramp;
            yaw_dot = yawrate;
        elseif t>=(tinit+tinter) && t<(tend-tinter)
            roll_dot = 0;
            pitch_dot = 0;
            yaw_dot = yawrate;
        elseif t>=(tend-tinter) && t<(tend-tinter+tramp)
            roll_dot = -rollrate*(t-(tend-tinter))/tramp;
            pitch_dot = -pitchrate*(t-(tend-tinter))/tramp;
            yaw_dot = yawrate;
        elseif t>=(tend-tinter+tramp) && t<(tend-tramp)
            roll_dot = -rollrate;
            pitch_dot = -pitchrate;
            yaw_dot = yawrate;
        elseif t>=(tend-tramp) && t<tend
            roll_dot = -rollrate*(tend-t)/tramp;
            pitch_dot = -pitchrate*(tend-t)/tramp;
            yaw_dot = yawrate*(tend-t)/tramp;            
        else
            roll_dot = 0;
            pitch_dot = 0;
            yaw_dot = 0;
        end
        
    case 3 % Transition
        
        tinit = varargin{2};
        tend  = varargin{3};
        tramp = varargin{4};
        tinter1 = varargin{5};
        tinter2 = varargin{6};
        
        if t>=tinit && t<(tinit+tramp)
            roll_dot = 0;
            pitch_dot = pitchrate*(t-tinit)/tramp;
            yaw_dot = 0;
        elseif t>=(tinit+tramp) && t<(tinit+tinter1-tramp)
            roll_dot = 0;
            pitch_dot = pitchrate;
            yaw_dot = 0;
        elseif t>=(tinit+tinter1-tramp) && t<(tinit+tinter1)
            roll_dot = 0;
            pitch_dot = pitchrate*((tinit+tinter1)-t)/tramp;
            yaw_dot = 0;
        elseif t>=(tinit+tinter1) && t<(tend-tinter2)
            roll_dot = 0;
            pitch_dot = 0;
            yaw_dot = 0;            
        elseif t>=(tend-tinter2) && t<(tend-tinter2+tramp)
            roll_dot = 0;
            pitch_dot = -pitchrate*(t-(tend-tinter2))/tramp;
            yaw_dot = 0;
        elseif t>=(tend-tinter2+tramp) && t<(tend-tramp)
            roll_dot = 0;
            pitch_dot = -pitchrate;
            yaw_dot = 0;
        elseif t>=(tend-tramp) && t<tend
            roll_dot = 0;
            pitch_dot = -pitchrate*(tend-t)/tramp;
            yaw_dot = 0;            
        else
            roll_dot = 0;
            pitch_dot = 0;
            yaw_dot = 0;
        end
        
    otherwise
        
        roll_dot = 0;
        pitch_dot = 0;
        yaw_dot = 0;
        
end

euler_dot = [
    roll_dot
    pitch_dot
    yaw_dot
    ];

end