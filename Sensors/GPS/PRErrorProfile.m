function PRerror = PRErrorProfile(time,interval,type,param)

T = time(2) - time(1);
N = numel(time);

PRerror = zeros(1,N);

tinit = interval(1);
tend = interval(2);

init = fix(tinit/T) + 1;
endi = fix(tend/T) + 1;

a = param.failure.a;
b = param.failure.b;

switch type
    case 'step'
        for i=init:endi
            PRerror(i) = b;
        end
    case 'ramp'
        for i=init:endi
            PRerror(i) = b/a*(i-init);
        end
    case 'sinusoidal'
        
end

end