function [Tr,stdTr] = TropoCorr(lat,H,el,D,param)
% DO-316 J.4

k1 = 77.604;   % K/mbar
k2 = 38200;    % K2/mbar
Rd = 287.054;  % J/(kg.K)
gm = 9.784;    % m/s
g = 9.8065;    % m/s
stdTVE = 0.12; % m

[P,T,e,beta,lambda] = zetaLatD(lat,D,param);

z_hyd = 1e-6*k1*Rd*P/gm;
z_wet = 1e-6*k2*Rd/(gm*(lambda+1)-beta*Rd)*e/T;

d_hyd = (1-beta*H/T)^(g/beta/Rd)*z_hyd;
d_wet = (1-beta*H/T)^(g*(lambda+1)/beta/Rd-1)*z_wet;

m_el = 1.001/sqrt(0.002001+sin(el)^2)*(1+0.015*max(0,4-el*180/pi));

Tr = -(d_hyd+d_wet)*m_el;
stdTr = stdTVE*m_el;

end

function [P,T,e,beta,lambda] = zetaLatD(lat,D,param)

Table0 = param.sensors.GPS.tropo_corr.Table0;
TableDelta = param.sensors.GPS.tropo_corr.TableDelta;

if lat>0
    Dmin = 28;
elseif lat<0
    Dmin = 211;
end

atmParam = NaN(4,1);
for i=1:5
    [zeta0,Dzeta] = linInterp(lat,Table0(:,i),TableDelta(:,i));
    atmParam(i) = zeta0-Dzeta*cos(2*pi*(D-Dmin)/365.25);
end

P      = atmParam(1);
T      = atmParam(2);
e      = atmParam(3);
beta   = atmParam(4);
lambda = atmParam(5);

end

function [zeta0,Dzeta] = linInterp(lat,z0,Dz)

if lat <= 15
    row0 = 1;
    row1 = 1;
    lat0 = 15;
    lat1 = 30;
elseif lat > 15 && lat <= 30
    row0 = 1;
    row1 = 2;
    lat0 = 15;
    lat1 = 30;
elseif lat > 30 && lat <= 45
    row0 = 2;
    row1 = 3;
    lat0 = 30;
    lat1 = 45;
elseif lat > 45 && lat <= 60
    row0 = 3;
    row1 = 4;
    lat0 = 45;
    lat1 = 60;
elseif lat > 60 && lat <= 75
    row0 = 4;
    row1 = 5;
    lat0 = 60;
    lat1 = 75;
elseif lat > 75
    row0 = 5;
    row1 = 5;
    lat0 = 75;
    lat1 = 75;
end

zeta0 = z0(row0)+(z0(row1)-z0(row0))*(lat-lat0)/(lat1-lat0);
Dzeta = Dz(row0)+(Dz(row1)-Dz(row0))*(lat-lat0)/(lat1-lat0);

end