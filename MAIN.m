close all
clear all
clc

deg2rad = pi/180;
rad2deg = 180/pi;
degph2radps=pi/180/3600;

ft2m = 0.3048;
m2ft = 1/ft2m;
mps2kn = 1.9438444924406;
kn2mps = 1/mps2kn;

% Parameters
mu = 3.986004418e14;      % [m3/s2]
OmegaE = 7.2921151467e-5; % [rad/s]
a = 6378137.0;           % [m]
b = 6356752.3;            % [m]
% f = 1/298.25;             % [] Earth flatness
g0 = 9.780327;            % [m/s2]
g_MEAN = 9.80665;         % [m/s2]

R0 = a;
f  = (a-b)/a;

c = 2.99792458e8;     % [m/s] speed of light
f_L1 = 1575.42e6;     % [Hz] L1 frequency
%OrbitData            % GNSS constellation

param.earth = struct('mu',mu,'OmegaE',OmegaE,'R0',R0,'f',f,'g0',g0,'c',c);

%% ============================= SETTINGS =================================
MNSsettings
% -------------------------------------------------------------------------

disp('### STARTING SIMULATION ###')
disp('Selected settings:')
fprintf('-> Flight profile: %s \n',FP)
fprintf('-> Kalman Filter: %s \n',filter)
fprintf('-> Filter corrections starts at: %3.1f s \n',correction_time)
fprintf('-> Total States Propagation starts at: %3.1f s \n\n',TSP_start)
if ChiSquareTestOff
    str = "OFF";
else
    str = "ON";
end
fprintf('-> Chi Square Test: %s \n', str)
if Underweighting
    str = "ON";
else
    str = "OFF";
end
fprintf('-> Innovation with underweighting: %s \n\n', str)

NsensorsList = numel(sensorsList);
for n=1:NsensorsList
    sensorName = sensorsList{n};
    eval([sensorName '_enabled = false;'])
    eval([sensorName '_DATAtype =' char(39) 'SIM' char(39) ';'])
end

Nsensors = numel(sensorSelected);
for n=1:Nsensors
    sensorName = sensorSelected{n};
    eval([sensorName '_enabled = true;'])
end

NFT = numel(FTList);
for n=1:NFT
    sensorName = FTList{n};
    eval([sensorName '_DATAtype =' char(39) 'FT' char(39) ';'])
end

% Filter selection
feedforward_enabled = false;
feedback_enabled = false;
calibration_enabled = false;
switch filter
    case 'LKF'
        feedforward_enabled = true;
    case 'EKF'
        feedforward_enabled = true;
        feedback_enabled = true;
        calibration_enabled = true;
    otherwise
        warning('No filter selected')
end

%% ================== LOAD DATA | SET SENSORS PARAMS ======================
SetSensorsParams

%% ===================== MEMORY PRE-ALLOCATION ============================

% ### TOTAL STATE ###
% X_INS = [lat long alt VN VE VD q0_bp q1_bp q2_bp q3_bp csi armX armY
% armZ]'
nXINS = 10+1+3;
% X_IMU = [nablaN nablaE nablaD epsilonN epsilonE epsilonD]';
nXIMU = 6;
% Y_INS = [PosX_ECEF PosY_ECEF PosZ_ECEF PosX_NED PosY_NED PosZ_NED 
% roll pitch yaw g_NED Asp_b_CG omega_bv_b]'
nY = 3+3+3+1+3+3;

nX = nXINS + nXIMU;
if GPSTC_enabled
    nX = nX + nGPS;
end
if ALT_enabled
    nX = nX + nALT;
end

% Kvert = [10 10 2]; % Gains for vertical channel
Kvert = [0 0 0];

% ### KALMAN FILTER ###
% x_NAV = [DRN DRE DRD DVN DVE DVD psiN psiE psiD]';
nxNAV = 9; 
% x_IMU = [DnablaN DnablaE DnablaD DepsilonN DepsilonE DepsilonD]';
nxIMU = 6; 

nKF = nxNAV + nXIMU;
Q_KF = Q_IMU;
P0_KF = P0_IMU;
if GPSTC_enabled
    nKF = nKF + nGPS;
    Q_KF = blkdiag(Q_KF,Q_GPS);
    P0_KF = blkdiag(P0_KF,P0_GPS);
end
if ALT_enabled
    nKF = nKF + nALT;
    Q_KF = blkdiag(Q_KF,Q_ALT);
    P0_KF = blkdiag(P0_KF,P0_ALT);
end

%  Error error state
xi_til = NaN(nKF,N,M);
x_til = NaN(nKF,N);
y_til = NaN(9,N);

% Time frames when update occured
kU = NaN(N,M);

% MC statistics
NEES = NaN(N,M);
NEES_TR = NaN(N,M);
NEES_NAV = NaN(N,M);
NEES_NI = NaN(N,M);
NEES_IMU = NaN(N,M);
NEES_GPS = NaN(N,M);
NIS = NaN(N,M);
DOF = NaN(N,M);

%% ======================== MONTE CARLO SIMULATION ========================

tic
if M==1
    seed = 52821;
    rng(seed);
else
    rng('shuffle');
end
disp(['Starting Monte Carlo simulation with ',num2str(M),' realization(s) ...'])
for i=1:M
    fprintf('===============> Realization #%d...\n',i)
    k = 1;

    % --------------------- MEMORY PRE-ALLOCATION -------------------------
    % Total states
    X = NaN(nX-2,N);%<-----------------------------------------------------
    Y = NaN(nY,N);
    % Estimated Total State
    X_hat = NaN(nX-2,N);%<-------------------------------------------------
    Y_hat = NaN(nY-6,N); % removes [Asp_b_CG omega_bv_b]'
    % True Error State
    x = NaN(nKF,N);
    y = NaN(9,N);
    % Estimated Error State
    x_hat = NaN(nKF,N);
	y_hat = NaN(9,N);
    % Estimated Covariance Matrix
    P = NaN(nKF,nKF,N);
    % Monitoring variable
    tracePnorm = NaN(N,1);
    
    % -------------------------- LOAD SENSORS -----------------------------
    
    % IMU
    if strcmp(IMU_DATAtype,'SIM')
        IMU_data = IMUi_data(:,:,i)';
        dAsp_b_m = IMU_data(1:3,:);
        dOmega_bi_b_m = IMU_data(4:6,:);
        IMUbias_GT = IMUibias_GT(:,1:3,i)';
        IMUdrift_GT = IMUibias_GT(:,4:6,i)';
    elseif strcmp(IMU_DATAtype,'FT')
        if strcmp(param.sensors.IMU.IMUoutput,'standard')
            Asp = IMU_data(1:3,:);
            Omega = IMU_data(4:6,:);
            IMUbias_GT = zeros(size(Asp));
            IMUdrift_GT = zeros(size(Omega));
        elseif strcmp(param.sensors.IMU.IMUoutput,'increments')
            dAsp_b_m = IMU_data(1:3,:)*T_IMU;
            dOmega_bi_b_m = IMU_data(4:6,:)*T_IMU;
            IMUbias_GT = zeros(size(dAsp_b_m));
            IMUdrift_GT = zeros(size(dOmega_bi_b_m));
        end
    end
    
    for n=1:Nsensors
        sensorName = sensorSelected{n};
        % Initialize nem measurement and flag variable
        eval(['new' sensorName ' = false;'])
        eval(['flag' sensorName ' = false;'])
        % Loading Monte Carlo Realizations
        if strcmp(eval([sensorName '_DATAtype']),'SIM')
            eval([sensorName '_data =' sensorName 'i_data{i};'])
            if strcmp(sensorName,'GPSTC')
                GPSTC_data = [GPSTC_data,GPSTC_aux];
            end
        end
        % Read sensors
        if eval([sensorName '_enabled'])
            eval(['m' sensorName '= 1;'])
            eval(['[Z_' sensorName ',new' sensorName ',flag' sensorName...
                '] = ' sensorName '_READ(' sensorName '_data,m' sensorName...
                ',param);'])
        end
    end
    
    % Bias state GT
    if GPSTC_enabled && strcmp(GPSTC_DATAtype,'SIM')
        GPSTCbias_GT = GPSTCibias_GT{i};
    elseif ~GPSTC_enabled
        GPSTCbias_GT = [];
    end
    if ALT_enabled && strcmp(ALT_DATAtype,'SIM')
        ALTbias_GT = ALTibias_GT{i};
    elseif ~ALT_enabled
        ALTbias_GT = [];
    end
    
    % ---------------------- INS INITIALIZATION ---------------------------
    INSinit

    % -------------------- KALMAN FILTER INITIALIZATION -------------------   
    % ### Initial Estimated Total State ###
    X_hat(:,k) = X0;
    Y_hat(:,k) = Y0(1:10);
    
    update = false;
    reset = false;
    lach_update = false;
    latch25 = true;
    
    if starWithFeedback
        feedback_on = true;
        feedback_latch = false;
        feedback_start = 1;
    else
        feedback_on = false;
        feedback_latch = true;
    end
    
    sw_param = [0 1];
    
    % ### Estimated Error State and Covariance Matrix###
    x0_hat = zeros(nKF,1);
    x0_hat(16:17) = [b0; 100];
    x_hat(:,k) = mvnrnd(x0_hat,P0_KF); % Initial estimated error vector
    P(:,:,k) = P0_KF;

    % EKF state matrix input
    U_EKF = zeros(6,1);
    
    y_hat(:,k) = zeros(9,1);

    % -------------------- KALMAN FILTER INITIAL UPDATE ------------------- 
    ChiSquareSum = 0;
    ChiDofSum = 0;
    for n=1:Nsensors
        sensorName = sensorSelected{n};
        if eval(['new' sensorName ' && flag' sensorName])
            % Sensor Measurement Model
            eval(['[z,H,R,aux] = ' sensorName '_MM(Z_' sensorName ',X(:,1),U_EKF,param);'])
            % EKF update
            [x_hat(:,k),P(:,:,k),ChiSquare,dof_z,update] = ...
                EKFupdate(k,x_hat(:,k),P(:,:,k),z,H,R,aux,param,feedback_on);
            eval(['new' sensorName ' = false;'])
            if ~lach_update && update
                lach_update = true;
                reset = true;
                kU(k,i) = k;
            end
            % ChiSquare store
            ChiSquareSum = ChiSquareSum + ChiSquare;
            ChiDofSum = ChiDofSum + dof_z;
        end
    end
    if lach_update
        NIS(k,i) = ChiSquareSum;
        DOF(k,i) = ChiDofSum;
    end

    % ------------------------- TRUE ERROR STATE --------------------------
    [x(1:9,k),y(:,k)] = TrueErrorState(X(:,k),X_GT(1,:)',Y_GT(1,:)',param);
    
    x(10:12,k) = X_hat(15:17,k) - IMUbias_GT(:,k);
    x(13:15,k) = X_hat(18:20,k) - IMUdrift_GT(:,k); 
    if GPSTC_enabled
        x(16:17,k) = GPSTCbias_GT(:,k);
    end
    if ALT_enabled
        x(end,k) = X_hat(end,k) - ALTbias_GT(k);
    end

    % ------------------------- ERROR ERROR-STATE -------------------------
    xi_til(:,k,i) = x_hat(:,k) - x(:,k);
    y_til(:,k) = y_hat(:,k) - y(:,k);
    
    NEES(k,i) = xi_til(:,k,i)'/P(:,:,k)*xi_til(:,k,i);
    NEES_TR(k,i) = sum(xi_til(:,k,i).*xi_til(:,k,i)./diag(P(:,:,k)));
    NEES_NAV(k,i) = xi_til(1:9,k,i)'/P(1:9,1:9,k)*xi_til(1:9,k,i);

    % ------------------------- FEEDFORWARD ---------------------------
    if ~feedback_on || lach_update
        [X_hat(:,k),Y_hat(:,k),y_hat(:,k)] = StateCorrection(X(:,k),x_hat(:,k),param);
    else
        X_hat(:,k) = X(:,k);
        Y_hat(:,k) = Y(1:10,k);
    end

    % Trace P
    tracePnorm(k) = sum(diag(P(:,:,k))./diag(P0_KF));%trace(P(:,:,k,i)/P0_KF);
    
    % --------------------- FEEDBACK AND RESET ------------------------
    if feedback_enabled && reset && feedback_on
        
        X(:,k) = X_hat(:,k);
        Y(1:10,k) = Y_hat(:,k);
        
%         x_hat(:,k) = zeros(nKF,1); <---------------------------------
        x_hat(1:15,k) = zeros(15,1);
        if ALT_enabled
            x_hat(end,k) = 0;
        end
        
        lach_update = false;
        reset = false;
    end
    
    %% ========================== TIME LOOP ===============================
    k0=2;
    for k=k0:N
        
        % ------------------------ SENSORS READING ------------------------
        for n=1:Nsensors
            sensorName = sensorSelected{n};
            % Read sensors
            if eval([sensorName '_enabled && mod(k-1,k' sensorName ')==0'])
                eval(['m' sensorName '= m' sensorName ' +1;'])
                eval(['[Z_' sensorName ',new' sensorName ',flag' sensorName...
                    '] = ' sensorName '_READ(' sensorName '_data,m' sensorName...
                    ',param);'])
            end
        end
        
        % ------------------------- TOTAL STATE ---------------------------
        if strcmp(param.sensors.IMU.IMUoutput,'increments')
            % 4 Samples Salichev Algorithm
            idx_IMU = (k-1)*4+1;
            U_IMU = [
                dAsp_b_m(:,idx_IMU-4:idx_IMU-1)
                dOmega_bi_b_m(:,idx_IMU-4:idx_IMU-1)
                ];
            % ### ONLINE CALIBRATION ###
            if calibration_enabled && feedback_on
                U_hat = U_IMU - X_hat(15:20,k-1)*T_IMU;
            else
                U_hat = U_IMU;
            end
            % altimeter = alt_GT(idx_GT);
            U = U_hat;
        elseif strcmp(IMUoutput,'standard')
            % RK4
            g_NED = Y(10,k-1);
            U_IMU = [
                Asp(:,k-1)*g_NED
                Omega(:,k-1)
                ];
            % ### ONLINE CALIBRATION ###
            if calibration_enabled && feedback_on
                U_hat = U_IMU - X_hat(15:20,k-1);
            else
                U_hat = U_IMU;
            end
            U = [U_hat; 0];
        end
            
        if time(k)>TSP_start %&& V_G > 0.01
            [X(:,k),Y(:,k)] = TotalStatePropagation(k,X(:,k-1),U,param); 
        else
            X(:,k) = X(:,k-1);
            Y(:,k) = Y(:,k-1);
        end
        
        % ------------------ KALMAN FILTER PREDICITION --------------------
        % Discrete Error-state Kalman Filter
        % ### Filter Prediction ###
        if strcmp(param.sensors.IMU.IMUoutput,'increments')
            U_EKF = sum(U,2)/T;
        elseif strcmp(param.sensors.IMU.IMUoutput,'standard')
            U_EKF = U;
        end
        [x_hat(:,k),P(:,:,k)] = EKFprediction(X_hat(:,k-1),U_EKF,...
            x_hat(:,k-1),P(:,:,k-1),Q_KF,param);
        
        % -------------------- KALMAN FILTER UPDATE -----------------------
        ChiSquareSum = 0;
        ChiDofSum = 0;
        lach_update = false;
        for n=1:Nsensors
            sensorName = sensorSelected{n};
            if eval(['new' sensorName ' && flag' sensorName])
                % Sensor Measurement Model
                eval(['[z,H,R,aux] = ' sensorName '_MM(Z_' sensorName...
                    ',X(:,k),U_EKF,param);']) %<-------------------------
                % EKF update
                [x_hat(:,k),P(:,:,k),ChiSquare,dof_z,update] = ...
                    EKFupdate(k,x_hat(:,k),P(:,:,k),z,H,R,aux,param,feedback_on);
                eval(['new' sensorName ' = false;'])
                if ~lach_update && update
                    lach_update = true;
                    reset = true;
                    kU(k,i) = k;
                end
                % ChiSquare store
                ChiSquareSum = ChiSquareSum + ChiSquare;
                ChiDofSum = ChiDofSum + dof_z;
            end
        end
        if lach_update
            NIS(k,i) = ChiSquareSum;
            DOF(k,i) = ChiDofSum;
        end
        
        % ----------------------- TRUE ERROR STATE ------------------------   
        % ### IDX ###
        idx_GT = fix((k-1)*ratioTdt+1);
        [x(1:9,k),y(:,k)] = TrueErrorState(X(:,k),X_GT(idx_GT,:)',Y_GT(idx_GT,:)',param);
        
        idx_IMU = fix((k-1)*ratioTTIMU+1);
        x(10:12,k) = X(15:17,k) - IMUbias_GT(:,idx_IMU);
        x(13:15,k) = X(18:20,k) - IMUdrift_GT(:,idx_IMU); 
        if GPSTC_enabled
            x(16:17,k) = GPSTCbias_GT(:,k);
        end
        if ALT_enabled
            x(end,k) = X(end,k) - ALTbias_GT(k);
        end
        
        % ----------------------- ERROR ERROR STATE -----------------------
        xi_til(:,k,i) = x_hat(:,k) - x(:,k);
        y_til(:,k) = y_hat(:,k) - y(:,k);
        
        NEES(k,i) = xi_til(:,k,i)'/P(:,:,k)*xi_til(:,k,i);
        NEES_TR(k,i) = sum(xi_til(:,k,i).*xi_til(:,k,i)./diag(P(:,:,k)));
        NEES_NAV(k,i) = xi_til(1:9,k,i)'/P(1:9,1:9,k)*xi_til(1:9,k,i);

        % ------------------------- FEEDFORWARD ---------------------------
        if ~feedback_on || lach_update
            [X_hat(:,k),Y_hat(:,k),y_hat(:,k)] = StateCorrection(X(:,k),x_hat(:,k),param);
        else
            X_hat(:,k) = X(:,k);
            Y_hat(:,k) = Y(1:10,k);
            y_hat(:,k) = y_hat(:,k-1);
        end
        
        % --------------------- FEEDBACK AND RESET ------------------------
        % Trace P
        tracePnorm(k) = sum(diag(P(:,:,k))./diag(P0_KF));
        difftrace = tracePnorm(k)-tracePnorm(k-1);
        if tracePnorm(k) < tracePnorm_threshold && difftrace<0 && ...
                feedback_latch && time(k)>=correction_time
            feedback_on = true;
            feedback_latch = false;
            feedback_start = k;
        end
        
        if feedback_enabled && reset && feedback_on
            X(:,k) = X_hat(:,k);
            Y(1:10,k) = Y_hat(:,k);
            x_hat(1:15,k) = zeros(15,1);
            if ALT_enabled
                x_hat(end,k) = 0;
            end
            lach_update = false;
            reset = false;
        end

        % ------------------------ SIM PROGRESS ---------------------------
        evolution = round(k/N*100);
        if evolution==25 && latch25
            fprintf('---------------> %d%%...\n',evolution)
            latch25 = false;
            latch50 = true;
        elseif evolution==50 && latch50
            fprintf('---------------> %d%%...\n',evolution)
            latch50 = false;
            latch75 = true;
        elseif evolution==75 && latch75
            fprintf('---------------> %d%%...\n',evolution)
            latch75 = false;
            latch100 = true;
        elseif evolution==100 && latch100
            fprintf('---------------> %d%%...\n',evolution)
            latch100 = false;
        end
    end
end
x_til(:,:) = xi_til(:,:,i);
toc
%% ============================== STATISTICS ==============================
if M>1
    tic
    disp('Running Monte Carlo Statistics')
    
    % Expected value
    E_x = sum(xi_til,3)/M;

    % Covariance
    error_x = NaN(nKF,N,M);
    Cov_x = zeros(nKF,nKF,N);
    for i=1:M
        error_x(:,:,i) = xi_til(:,:,i) - E_x;
    end
    for k=1:N
        for i=1:M
            Cov_x(:,:,k) = Cov_x(:,:,k) + error_x(:,k,i)*error_x(:,k,i)';
        end
    end
    Cov_x = Cov_x/(M-1);
    
    % Consistency checks
    ANEES = sum(NEES,2);
    ANEES_TR = sum(NEES_TR,2)/M;
    ANEES_NI = sum(NEES_NAV,2)/M;
    ANEES_NAV = sum(NEES_NAV,2)/M;
    ANEES_IMU = sum(NEES_IMU,2)/M;
    ANEES_GPS = sum(NEES_GPS,2)/M;
    ANIS = sum(NIS,2)/M;
    NIS_DOF = 0;
    toc
else
    E_x = 0;
    error_x = 0;
    Cov_x = 0;
    ANEES = 0;
    ANEES_TR = 0;
    ANEES_NI = 0;
    ANEES_NAV = 0;
    ANEES_IMU = 0;
    ANEES_GPS = 0;
    ANIS = NIS;
end
%% =========================== DATA SAVING ================================

if saveData
    dateStamp = char(datetime('now','Format','MMdd_HHmm'));%'yyyyMMdd_HHmmss'
    Ns = numel(sensorSelected);
    cat_aux = '';
    for n=1:Ns
        cat_aux = strcat(cat_aux,'_',sensorSelected(n));
    end
    eval(['save SIM_' dateStamp '_' FP '_' IMUtype int2str(IMUmodel)...
        cat_aux{1} '_M' int2str(M) '.mat ' ...
        'time X Y X_hat Y_hat x y x_hat y_hat '...
        'P x_til y_til kU E_x Cov_x '...
        'ANEES ANEES_NI ANEES_NAV ANEES_IMU ANEES_GPS ANEES_TR ANIS ' ...
        'IMUbias_GT IMUdrift_GT GPSTCbias_GT ALTbias_GT '...
        'DOF tracePnorm sensorSelected FTList param -v7.3'])
end

    
    