function [x_hat,P,ChiSquare,dof_z,update] = EKFupdate(k,x_hat,P,z,H,R,aux,param,feedback_on)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

chi2per          = param.config.chi2per;
ChiSquareTestOff = param.config.ChiSquareTestOff;
ChiSquareGain    = param.config.ChiSquareGain;
UW_enabled       = param.config.Underweighting;

T = param.sim.T;

sensorName = aux{1};
UB = aux{2};

update = false;
if strcmp(sensorName,'GPSTC') || strcmp(sensorName,'LPS')
    m = size(H,1);
    chiSquare = zeros(m,1);
    for j=1:m      
        W = H(j,:)*P*H(j,:)';
        % beta
        if UW_enabled
            Pstar = diag(H(j,:)~=0)*P;
            beta = (UB(j)/2)*(trace(Pstar))^2/trace(W);
        else
            beta = 0;
        end
        % Innovation w/ underweghting
        S = (1+beta)*W + R(j,j);
        % S = H(j,:)*P*H(j,:)'+R(j,j);
        z_til = z(j) - H(j,:)*x_hat;
        % ChiSquare (Inovation) Test
        chiSquare(j) = z_til*z_til/S;
        if chiSquare(j) < chi2inv(chi2per,1)*ChiSquareGain || ChiSquareTestOff || ~feedback_on
            % Kalman gain
            K = P*H(j,:)'/S;
            % Update
            x_hat = x_hat + K*z_til;
            % P = P - K*(H(j,:)*P);
            nKF = numel(x_hat);
            P = (eye(nKF)-K*H(j,:))*P*(eye(nKF)-K*H(j,:))'+K*R(j,j)*K';
            update = true;
        else
            fprintf('Outlier from %s at t = %6.2f s\n',sensorName,(k-1)*T)
            if strcmp(sensorName,'GPSTC')
                SVID = [aux{3} aux{3}];
                fprintf('Faulty SVID: %d\n',SVID(j))
            end
        end
    end
    chol(P); % checking if P is positive definite
    P = (P+P')/2; % Forcing P to be symmetric
    ChiSquare = sum(chiSquare);
    dof_z = m;
elseif strcmp(sensorName,'CAMLMS') || strcmp(sensorName,'CAMLMF') 
    m = numel(z)/2;
    chiSquare = zeros(m,1);
    for j=1:m
        S = H(2*j-1:2*j,:)*P*H(2*j-1:2*j,:)'+R(2*j-1:2*j,2*j-1:2*j);
        z_til = z(2*j-1:2*j) - H(2*j-1:2*j,:)*x_hat;
        % ChiSquare (Inovation) Test
        chiSquare(j) = z_til'/S*z_til;
        if chiSquare(j) < chi2inv(chi2per,2)*ChiSquareGain || ChiSquareTestOff || ~feedback_on
            % Kalman gain
            K = P*H(2*j-1:2*j,:)'/S;
            % Update
            x_hat = x_hat + K*z_til;
            % P = P - K*(H(2*j-1:2*j,:)*P);
            nKF = numel(x_hat);
            P = (eye(nKF)-K*H(2*j-1:2*j,:))*P*(eye(nKF)-K*H(2*j-1:2*j,:))'+K*R(2*j-1:2*j,2*j-1:2*j)*K';
            update = true;
        else
            fprintf('Outlier from %s at t = %6.2f s\n',sensorName,(k-1)*T)
            LM = aux{3};
            fprintf('Faulty LM: %d/%d\n',LM,j)
        end
    end
    chol(P); % checking if P is positive definite
    P = (P+P')/2; % Forcing P to be symmetric
    ChiSquare = sum(chiSquare);
    dof_z = m;
else
    W = H*P*H';
    % beta
    if UW_enabled
        Pstar = diag(H(1,:)~=0)*P;
        beta = (UB/2)*(trace(Pstar))^2/trace(W);
    else
        beta = 0;
    end
    
    % Innovation w/ underweghting
    S = (1+beta)*W + R;
    z_til = z - H*x_hat;
    % ChiSquare (Inovation) Test
    ChiSquare = z_til'/S*z_til;
    
    dof_z = numel(z);
    if ChiSquare < chi2inv(chi2per,dof_z)*ChiSquareGain || ChiSquareTestOff || ~feedback_on
        % Kalman gain
        K = P*H'/S;
        % Update
        x_hat = x_hat + K*z_til;
        nKF = numel(x_hat);
        % P = (eye(nKF)-K*H)*P;
        P = (eye(nKF)-K*H)*P*(eye(nKF)-K*H)'+K*(beta*H*P*H'+R)*K'; % Joseph's formula w/ underweghting
        chol(P); % checking if P is positive definite
        P = (P+P')/2; % Forcing P to be symmetric
        update = true;
    else
        update = false;
        fprintf('Outlier from %s at t = %6.2f s\n',sensorName,(k-1)*T)
    end
end

end