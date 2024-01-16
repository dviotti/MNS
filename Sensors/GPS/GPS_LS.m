function [R_rec,b] = GPS_LS(R_SV,rho_m,varargin)

[p,q] = size(rho_m);
if p==1
    rho_m = rho_m';
    p = q;
elseif ~(p>=1 && q==1)
    error('Second argument must be a nx1 vector')
end

[~,s]= size(R_SV);
if p~=s
    error('Number of visible SVs is different from pseudoresnge measurements')
elseif s<4
    warning('Number of visible SVs is less than 4. No solution can be found')
    R_rec = NaN(3,1);
    b = NaN;
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

x = zeros(4,1);
Dx = 10*ones(4,1);
eps = 1;%1e-3;

i=1;
while norm(Dx)>eps

    R_rec = x(1:3);
    b = x(4);
    
    [LOS,dist] = GPSlinear(R_SV,R_rec);
    H = [-LOS ones(p,1)];
     
    rho_hat = dist + b;
    Drho = rho_m - rho_hat;

    Hstar = (H'*W*H);
    if det(Hstar) <1e-20
        warning('Hstar ill-conditioned. GPS solution not found')
        R_rec = NaN(3,1);
        b = NaN;
        return
    end
        
    A = eye(4)/Hstar;
    Dx = A*H'*W*Drho;
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
        R_rec = NaN(1:3);
        b = NaN;
        return
    end
    i = i+1;
end

R_rec = x(1:3);
b = x(4);

end