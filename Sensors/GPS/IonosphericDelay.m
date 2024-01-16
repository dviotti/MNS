%% Ionospheric delay
%
% [T_IONO, D_IONO, STD_IONO, VAR_IONO] = IonosphericDelay(GPSTIME,LATITUDE,LONGITUDE,AZIMUT,ELEVATION,ALPHA,BETA)
%  returns:
%   - Ionospheric delay(s) T_IONO [s],
%   - Ionospheric delay(s) D_IONO [m],
%   - Ionospheric standard deviation(s) STD_IONO [m],
%   - Ionospheric variance(s) VAR_IONO [m²],
%  from:
%   - GPS time of week GPSTIME [s],
%   - Position of the receiver, in terms of LATITUDE and LONGITUDE [°],
%   - Satellite LOS, in terms of AZIMUT and ELEVATION [°],
%   - 4-alpha and 4-beta coefficients, of Klobuchar ionospheric model, in terms of vectors ALPHA and BETA.
%
% Sources:
%  - Ionospheric delay    : IS-GPS-200 20.3.3.5.2.2
%  - Ionospheric variance : DO-229 J.2.3
%
% Example:
%  GPStime    = 0;
%  Latitude   = 0;
%  Longitude  = 0;
%  alpha      = [+2.887e-08, +2.980e-08, -1.192e-07,       0];
%  beta       = [   +153600,    -196608,     -65536, +393216];
%  Azimuts    = 360*rand(1,5);
%  Elevations =  90*rand(1,5);
%  [T_iono, D_iono, Std_iono, Var_iono] = IonosphericDelay(GPStime,Latitude,Longitude,Azimuts,Elevations,alpha,beta);
%  fprintf(1,'\n %s | %s | %s \n','LOS','Delay [m]','STD [m]');
%  fprintf(1,'%s\n',repmat('-',1,27));
%  arrayfun(@(n)fprintf(1,' %3u | %9.3f | %7.3f\n',n,D_iono(n),sqrt(Var_iono(n))),1:numel(D_iono));

function [T_iono, D_iono, Std_iono, Var_iono] = IonosphericDelay(GPStime,Latitude,Longitude,Azimut,Elevation,alpha,beta)

% Control
if ~isscalar(GPStime)   || ~isa(GPStime,  'double'), error('IonosphericDelay: Invalid GPS time.');                  end
if ~isscalar(Latitude)  || ~isa(Latitude, 'double'), error('IonosphericDelay: Invalid latitude.');                  end
if ~isscalar(Longitude) || ~isa(Longitude,'double'), error('IonosphericDelay: Invalid longitude.');                 end
if ~all(eq(size(Azimut),size(Elevation))),           error('IonosphericDelay: Inconsistent azimut and elevation.'); end
if ne(numel(alpha),4),                               error('IonosphericDelay: Invalid vector ''alpha''.');          end
if ne(numel(beta),4),                                error('IonosphericDelay: Invalid vector ''beta''.');           end

% Constants
c  = 299792458.0;	% Light celerity        [m/s]
Re = 6378136;       % Earth radius          [m]
hI = 350000.0;      % Ionosphere altitude   [m]

% Conversion functions
deg2rad = @(a)pi/180*a;
rad2sc  = @(a)a/pi;
sc2rad  = @(a)pi*a;

% Possible transposition of Klobuchar coefficients for vectorial computation
if iscolumn(alpha), alpha = alpha'; end
if iscolumn(beta),  beta  = beta';  end

% Inputs
phi_u    = deg2rad(Latitude);       % User latitude         [rad]
lambda_u = deg2rad(Longitude);      % User longitude        [rad]
A        = deg2rad(Azimut);     	% Satellite azimut      [rad]
E        = deg2rad(Elevation);  	% Satellite elevation   [rad]

% Psi
Psi = 0.0137./(rad2sc(E)+0.11)-0.022;

% Pierce point latitude
phi_i = rad2sc(phi_u) + Psi.*cos(A);
phi_i = min(phi_i,+0.416);
phi_i = max(phi_i,-0.416);

% Pierce point longitude
lambda_i = rad2sc(lambda_u) + Psi .* sin(A)./cos(sc2rad(phi_i));

% Time
t = 4.32e4*lambda_i+GPStime;
t = t + 86400*(lt(t,0)-ge(t,86400));

% Geomagnetic latitude [sc]
phi_m = phi_i + 0.064*cos(sc2rad(lambda_i-1.617));

% Obliquity factor
F = 1.0 + 16.0 * (0.53-rad2sc(E)).^3;

% Period [s]
if iscolumn(phi_m)
    PER = (beta * phi_m.^(0:3)')';
else
    PER = beta * phi_m'.^(0:3)';
end
PER = max(PER,72000);

% Phase [rad]
x = 2*pi*(t-50400)./PER;

% Delay [s]
if iscolumn(phi_m)
    AMP = (alpha * phi_m.^(0:3)')';
else
    AMP = alpha * phi_m'.^(0:3)';
end
AMP = max(AMP,0);

% Ionospheric correction [s] (IS-GPS-200 20.3.3.5.2.2)
% (if the user is operating on the L2 frequency, the correction term must
% be multiplied by gamma)
T_iono = F * 5.0e-9;
T_iono = T_iono + lt(abs(x),1.57) .* F .* AMP .* (1-x.^2/2 + x.^4/24);

% Vertical deviation [m] (DO-229 J.2.3)
tau_vert = 6*ones(size(Azimut));
tau_vert(le(180*abs(phi_m),55)) = 4.5;
tau_vert(le(180*abs(phi_m),20)) = 9;

% Obliquity factor (DO-229 A.4.4.10.4)
Fpp = (1-(Re/(Re+hI)*cos(E)).^2).^-0.5;

% Variance of the ionospheric correction [m²] (DO-229 J.2.3)
Var_iono = max((c*T_iono/5).^2,(Fpp.*tau_vert).^2);

% Additional data
D_iono   = c*T_iono;
Std_iono = sqrt(Var_iono);

end
