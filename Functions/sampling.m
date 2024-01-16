function Y = sampling(U,dt,T)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

k = T/dt;

[n,m] = size(U);
N = fix(n/k)+1;
Y = zeros(N,m);
for i=1:N
    Y(i,:) = U(fix((i-1)*k)+1,:);
end

end

