function [R_SV_eR,R_rec,b] = GPSsolution(R_SV_eT,rho_m,URA_m,param)

% R_SV_etT: SV position in the ECEF at transmition time
% R_SV_etR: SV position in the ECEF at reception time

c = param.earth.c;
OmegaE = param.earth.OmegaE;

visibleSVs = length(rho_m);

R_SV_eR = zeros(3,visibleSVs);

% Correcting transmission ECEF frame
Dt = zeros(visibleSVs,1);
for j=1:visibleSVs
    R_SV_eR(:,j) = DCM(3,OmegaE*Dt(j))*R_SV_eT(:,j);
end
[R_rec,b] = GPS_LS(R_SV_eR,rho_m,URA_m);

Dt_prev = Dt;
range_m = rho_m - b*ones(size(rho_m));
Dt = range_m./c;

i = 1;
while abs(Dt-Dt_prev)>1e-6 % remove the while loop? 1 iteration is enough

    for j=1:visibleSVs
        R_SV_eR(:,j) = DCM(3,OmegaE*Dt(j))*R_SV_eT(:,j);
    end
    [R_rec,b] = GPS_LS(R_SV_eR,rho_m);

    Dt_prev = Dt;
    range_m = rho_m - b*ones(size(rho_m));
    Dt = range_m./c;

    if i>10
        warning('GPS solution not found')
        break
    end
    i = i + 1;
    % disp('oi')
end

end