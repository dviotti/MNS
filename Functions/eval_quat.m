function q = eval_quat(method,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if method==1
    phi = varargin{1};
    u = varargin{2};
    q = [
        cos(phi/2)
        u*sin(phi/2)
        ];
elseif  method==2
    angle_vec = varargin{1};
    phi = norm(angle_vec);
    if abs(phi)<1e-6
        u = [0 0 0]';
    else
        u = angle_vec/phi;
    end
    q = [
        cos(phi/2)
        u*sin(phi/2)
        ];
else
    error('Choose a method to evaluate the quaternion')
end

q = q/norm(q);

end

