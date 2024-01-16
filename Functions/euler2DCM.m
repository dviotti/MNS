function D = euler2DCM(euler,sequence)

[n,m] = size(euler);
if ~(n==3 && m==1 || n==1 && m==3)
    error('The first input must be 3x1 vector or 1x3 list')
end

if ~ischar(sequence)
    error('The second input must be a string')
end

switch sequence
    case 'ZYX'
        D = DCM(1,euler(1))*DCM(2,euler(2))*DCM(3,euler(3));
    case 'ZYZ'
        D = DCM(3,euler(1))*DCM(2,euler(2))*DCM(3,euler(3));
    otherwise
        error('No existent sequence selected')
end

end