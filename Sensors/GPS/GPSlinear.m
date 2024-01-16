function [LOS,range] = GPSlinear(R_SV,R_rec)

n = size(R_SV,2);

range = zeros(n,1); % Distance betweent Rx and SVi
LOS = zeros(n,3);  % Line-of-sigth vector from SVi to Rx
for i=1:n
    range(i) = sqrt((R_SV(1,i)-R_rec(1))^2+(R_SV(2,i)-R_rec(2))^2+(R_SV(3,i)-R_rec(3))^2);
    for j=1:3   
        LOS(i,j) = (R_SV(j,i)-R_rec(j))/range(i);
    end
end

end