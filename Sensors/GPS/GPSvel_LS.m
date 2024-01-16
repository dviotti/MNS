function [V_rec_NED,b_dot] = GPSvel_LS(R_SV,LLA_rec,V_SV,vel_m,varargin)

[p,q] = size(vel_m);
if p==1
    vel_m = vel_m';
    p = q;
elseif ~(p>=1 && q==1)
    error('Second argument must be a nx1 vector')
end

[~,s]= size(V_SV);
if p~=s
    error('Number of visible SVs is different from delta range measurements')
elseif s<4
    warning('Number of visible SVs is less than 4. No solution can be found')
    V_rec_NED = NaN(3,1);
    b_dot = NaN;
    return
end

if ~isempty(varargin)
    weight = varargin{1};
    if numel(weight)~=p
        error('Number of weight parameters is different from pseudoresnge measurements')
    end
    % W = diag(1./(weight.^2));
    W = diag(1./weight);
else
    W = eye(p);
end

lat = LLA_rec(1);
long = LLA_rec(2);
alt = LLA_rec(3);

R_rec = LLA2ECEF(lat,long,alt);
[LOS,~] = GPSlinear(R_SV,R_rec);
H = [-LOS ones(p,1)];

x = zeros(4,1);
Dx = 10*ones(4,1);
eps = 1;%1e-3;

i=1;
while norm(Dx)>eps

    V_rec = x(1:3);
    b_dot = x(4);
     
    vel_hat = diag(LOS*(V_SV - repmat(V_rec,1,s))) + b_dot;
    Dvel = vel_m - vel_hat;

    Hstar = (H'*W*H);
    if det(Hstar) <1e-20
        warning('Hstar ill-conditioned. GPS solution not found')
        V_rec_NED = NaN(3,1);
        b_dot = NaN;
        return
    end
        
    A = eye(4)/Hstar;
    Dx = A*H'*W*Dvel;
    x = x + Dx;
    % x'
    
    % R = chol(H'*H);
    % A = R\(R'\eye(4));
    % Dx = A*H'*Drho;
    % x = x + Dx;
    
%     Dx = LSwithQR(H,Drho);
%     x = x + Dx; 
    
    if i>15
        warning('GPS solution not found')
        V_rec_NED = NaN(1:3);
        b_dot = NaN;
        return
    end
    i = i+1;
end

V_rec = x(1:3);
b_dot = x(4);

D_ce = DCM(2,-(lat+pi/2))*DCM(3,long);
V_rec_NED = D_ce*V_rec;

end