function S = skew(X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%[n,m] = size(X);
eps = 10;

if iscolumn(X) %n==1 && m==3 || n==3 && m==1
    S = [0 -X(3) X(2)
        X(3) 0 -X(1)
        -X(2) X(1) 0];
%elseif n==3 && m==3
else
    isSkew = abs(norm(X+X')/2)<eps;
    if isSkew
        S1 = (X(3,2)-X(2,3))/2;
        S2 = (X(1,3)-X(3,1))/2;
        S3 = (X(2,1)-X(1,2))/2;
        S = [S1 S2 S3]';
    else
        error('Input is not a skew-symmetric matrix')
    end
%     S1 = (X(3,2)-X(2,3))/2;
%     S2 = (X(1,3)-X(3,1))/2;
%     S3 = (X(2,1)-X(1,2))/2;
%     S = [S1 S2 S3]';
%     if ~isSkew
%         warning('Input is not a skew-symmetric matrix')
%     end
% else
%     error('Input must be a 3D vector or a 3x3 matrix')
end

end

