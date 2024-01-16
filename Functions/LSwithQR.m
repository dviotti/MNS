function x = LSwithQR(A,b)

[m,n] = size(A);
if m<n
    error('The number of rows must be greater or equal than the number of columns')
end

[Q,R]=qr(A);

R1 = R(1:n,:);
bt = Q'*b;
c = bt(1:n);

x = R1\c;

end