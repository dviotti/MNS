function DCM = quat2DCM(method,varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% q = [1 0 0 0]';

if method==1
    q = varargin{1};
    [m,n] = size(q);
    if m==1 && n==4
        q = q';
    elseif ~(m==4 && n==1)
        error('Second argument must be a 4x1 quaternion')
    end
    q = q/norm(q);
elseif method==2
    q = eval_quat(2,varargin{1});
else
   error('Choose a method to evaluate the quaternion') 
end

qw = q(1);
q_vec = q(2:4);

% Right-hand rotation
% DCM = (qw^2-q_vec'*q_vec)*eye(3)+2*q_vec*(q_vec')+2*qw*skew(q_vec); 
% q_vec'*q_vec = 1 - qw^2;
% DCM = (2*qw^2-1)*eye(3)+2*q_vec*(q_vec')+2*qw*skew(q_vec);

% % Left-hand rotation
DCM = (qw^2-q_vec'*q_vec)*eye(3)+2*q_vec*(q_vec')-2*qw*skew(q_vec); 
% % q_vec'*q_vec = 1 - qw^2;
% % C = (2*qw^2-1)*eye(3)+2*q_vec*(q_vec')-2*qw*skew(q_vec);

end

