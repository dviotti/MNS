function VNED_dot = VNEDdotProfile(t,profile,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ANED  = varargin{1};
AN = ANED(1);
AE = ANED(2);
AD = ANED(3);

switch profile
    
    case 1 % Velocity change
        
        tinit = varargin{2};
        tend  = varargin{3};
        tramp = varargin{4};
        
        if t>=tinit && t<(tinit+tramp)
            VN_dot = AN*(t-tinit)/tramp;
            VE_dot = AE*(t-tinit)/tramp;
            VD_dot = AD*(t-tinit)/tramp;
        elseif t>=(tinit+tramp) && t<(tend-tramp)
            VN_dot = AN;
            VE_dot = AE;
            VD_dot = AD;
        elseif t>=(tend-tramp) && t<tend
            VN_dot = AN*(tend-t)/tramp;
            VE_dot = AE*(tend-t)/tramp;
            VD_dot = AD*(tend-t)/tramp;
        else
            VN_dot = 0;
            VE_dot = 0;
            VD_dot = 0;
        end
        
    case 2 % Position change (Takeoff and Landing)
        
        tinit = varargin{2};
        tend  = varargin{3};
        tramp = varargin{4};
        tinter = varargin{5};
        
        if t>=tinit && t<(tinit+tramp)
            VN_dot = AN*(t-tinit)/tramp;
            VE_dot = AE*(t-tinit)/tramp;
            VD_dot = AD*(t-tinit)/tramp;
        elseif t>=(tinit+tramp) && t<(tinit+tinter-tramp)
            VN_dot = AN;
            VE_dot = AE;
            VD_dot = AD;
        elseif t>=(tinit+tinter-tramp) && t<(tinit+tinter)
            VN_dot = AN*((tinit+tinter)-t)/tramp;
            VE_dot = AE*((tinit+tinter)-t)/tramp;
            VD_dot = AD*((tinit+tinter)-t)/tramp;
        elseif t>=(tend-tinter) && t<(tend-tinter+tramp)
            VN_dot = -AN*(t-(tend-tinter))/tramp;
            VE_dot = -AE*(t-(tend-tinter))/tramp;
            VD_dot = -AD*(t-(tend-tinter))/tramp;
        elseif t>=(tend-tinter+tramp) && t<(tend-tramp)
            VN_dot = -AN;
            VE_dot = -AE;
            VD_dot = -AD;
        elseif t>=(tend-tramp) && t<tend
            VN_dot = -AN*(tend-t)/tramp; 
            VE_dot = -AE*(tend-t)/tramp; 
            VD_dot = -AD*(tend-t)/tramp;            
        else
            VN_dot = 0;
            VE_dot = 0;
            VD_dot = 0;
        end   
        
    case 3 % N/S->E/W Coordinated Turn
        
        tinit = varargin{2};
        tend  = varargin{3};
        tramp = varargin{4};
        phase = varargin{5};
        
        AH = sqrt(AN^2+AE^2);
        
        if t>=tinit && t<(tinit+tramp)
            VN_dot = 0;
            VE_dot = sign(AE)*AH*cos(phase)*(t-tinit)/tramp;
            VD_dot = AD;
        elseif t>=(tinit+tramp) && t<(tend-tramp)
            VN_dot = sign(AN)*AH*sin((pi/2-phase)*(t-(tinit+tramp))/((tend-tramp)-(tinit+tramp))+phase);
            VE_dot = sign(AE)*AH*cos((pi/2-phase)*(t-(tinit+tramp))/((tend-tramp)-(tinit+tramp))+phase);
            VD_dot = AD;
        elseif t>=(tend-tramp) && t<tend
            VN_dot = sign(AN)*AH*sin(pi/2-phase)*(tend-t)/tramp;
            VE_dot = 0;
            VD_dot = AD;
        else
            VN_dot = 0;
            VE_dot = 0;
            VD_dot = 0;
        end
        
    case 4 % E/W->N/S Coordinated Turn
        
        tinit = varargin{2};
        tend  = varargin{3};
        tramp = varargin{4};
        phase = varargin{5};
        
        AH = sqrt(AN^2+AE^2);
        
        if t>=tinit && t<(tinit+tramp)
            VN_dot = sign(AN)*AH*cos(phase)*(t-tinit)/tramp;
            VE_dot = sign(AE)*AH*sin(phase)*(t-tinit)/tramp;
            VD_dot = AD;
        elseif t>=(tinit+tramp) && t<(tend-tramp)
            VN_dot = sign(AN)*AH*cos((pi/2-phase)*(t-(tinit+tramp))/((tend-tramp)-(tinit+tramp))+phase);
            VE_dot = sign(AE)*AH*sin((pi/2-phase)*(t-(tinit+tramp))/((tend-tramp)-(tinit+tramp))+phase);
            VD_dot = AD;
        elseif t>=(tend-tramp) && t<tend
            VN_dot = 0;
            VE_dot = sign(AE)*AH*(tend-t)/tramp;
            VD_dot = AD;
        else
            VN_dot = 0;
            VE_dot = 0;
            VD_dot = 0;
        end
        
    case 5 % Transition / Level change
        
        tinit = varargin{2};
        tend  = varargin{3};
        tramp = varargin{4};
        tinter = varargin{5};
        
        if t>=tinit && t<(tinit+tramp)
            VN_dot = AN*(t-tinit)/tramp;
            VE_dot = AE*(t-tinit)/tramp;
            VD_dot = AD*(t-tinit)/tramp;
        elseif t>=(tinit+tramp) && t<(tinit+tinter-tramp)
            VN_dot = AN;
            VE_dot = AE;
            VD_dot = AD;
        elseif t>=(tinit+tinter-tramp) && t<(tinit+tinter)
            VN_dot = AN;
            VE_dot = AE;
            VD_dot = AD*((tinit+tinter)-t)/tramp;
        elseif t>=(tinit+tinter) && t<(tend-tinter)
            VN_dot = AN;
            VE_dot = AE;
            VD_dot = 0;            
        elseif t>=(tend-tinter) && t<(tend-tinter+tramp)
            VN_dot = AN;
            VE_dot = AE;
            VD_dot = -AD*(t-(tend-tinter))/tramp;
        elseif t>=(tend-tinter+tramp) && t<(tend-tramp)
            VN_dot = AN;
            VE_dot = AE;
            VD_dot = -AD;
        elseif t>=(tend-tramp) && t<tend
            VN_dot = AN*(tend-t)/tramp;
            VE_dot = AE*(tend-t)/tramp;
            VD_dot = -AD*(tend-t)/tramp;            
        else
            VN_dot = 0;
            VE_dot = 0;
            VD_dot = 0;
        end
        
    otherwise
        
        VN_dot = 0;
        VE_dot = 0;
        VD_dot = 0;
        
end

VNED_dot = [
    VN_dot
    VE_dot
    VD_dot
    ];

end

