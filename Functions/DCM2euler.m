function euler = DCM2euler(D,sequence)

[n,m] = size(D);
if n~=3 && m~= 3
    error('The first input must be 3x3 matrix')
end

if ~ischar(sequence)
    error('The second input must be a string')
end

switch sequence
    case 'ZYX'
        phi   = atan2(D(2,3),D(3,3));
        theta = -asin(D(1,3));
        psi   = atan2(D(1,2),D(1,1));
    case 'ZYZ'
        phi   = atan2(D(2,3),-D(1,3));
        theta = acos(D(3,3));
        psi   = atan2(D(3,2),D(3,1));
    otherwise
        error('No existent sequence selected')
end

euler = [
    phi
    theta
    psi
    ];

end