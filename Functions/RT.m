function sigOUT = RT(sigIN,Tin,Tout,varargin)

if ~isempty(varargin)
    type = varargin{1};
else
    type = 'zoh';
end

[n,m] = size(sigIN);


if Tin>Tout
    k = fix(Tin/Tout);
    N = (n-1)*k+1;
    sigOUT = NaN(N,m);
    if strcmp(type,'zoh')
        for i=1:n-1
            idx1 = (i-1)*k+1;
            idx2 = i*k;
            sigOUT(idx1:idx2,:) = repmat(sigIN(i,:),k,1);
        end
        sigOUT(end,:) = sigIN(end,:);
    elseif strcmp(type,'foh')
        for i=1:n-1
            idx = (i-1)*k+1;
            sigOUT(idx,:) = sigIN(i,:);
            step = (sigIN(i+1,:)-sigIN(i,:))/k;
            for j=1:k
                sigOUT(idx+j,:) = sigOUT(idx+j-1,:)+step;
            end
        end
%         sigOUT(end-k:end,:) = repmat(sigIN(end,:),k,1);
    end
elseif Tin<Tout
    k = fix(Tout/Tin);
    sampling = 1:k:n;
    sigOUT = sigIN(sampling,:);
else
    sigOUT = sigIN;
end

end