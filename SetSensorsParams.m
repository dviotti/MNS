%% ================== LOAD DATA | SET SENSORS PARAMS ======================

% ### IMU ###
disp('Loading IMU data ...')
if strcmp(FP,'FTE2')
    load Sensors/IMU/IMU_data_FTE2.mat
else
    eval(['load Sensors/IMU/IMU_data_' FP '_' IMUtype '_' num2str(IMUmodel) '.mat'])
end
if strcmp(IMU_DATAtype,'SIM')
    param.sensors.IMU = paramIMU;
    param.AC.IMUarm = paramIMU.IMUarm;
    clear paramIMU
    % ### IMU OUTPUT ###
    param.sensors.IMU.IMUoutput = 'increments';
elseif strcmp(IMU_DATAtype,'FT')
    param.AC = paramAC;
    clear paramAC
    T_IMU = IMU_time(2) - IMU_time(1);

    switch param.AC.IMUtype
        case 'AHRS'
            IMUbias  	= 5*1e-3*g0;       % [m/s^2] - Accelerometer bias
            IMUdrift    = 180*degph2radps/10; % [rad/s] - RateGyro drift
            Accstd      = 1.5*1e-3*g0;     % [m/s^2] - Accelerometer standard deviation
            RateGyrostd = 27*degph2radps;  % [rad/s] - RateGyro standard deviation
        case 'IRS'
            IMUbias  	= 2.5*1e-3*g0;     % [m/s^2] - Accelerometer bias
            IMUdrift    = 36/10*degph2radps;  % [rad/s] - RateGyro drift
            Accstd      = 5/3*1e-3*g0;     % [m/s^2] - Accelerometer standard deviation
            RateGyrostd = 12*degph2radps;  % [rad/s] - RateGyro standard deviation
    end
    
    % ### IMU MODEL ###
    IMUmodel = 1; %1-constant bias; 3- Bias Random Walk; 3- Gauss-Markov; 4- GM no noise;
    % Constant Bias:    IMUtau = 1e12;  biasIns_on = 0;
    % Bias Random Walk: IMUtau = 1e12;  biasIns_on = 1;
    % Gauss-Markov:     IMUtau = 150;    biasIns_on = 1;
    IMUtau = 1e12;% [s] - Time constant
    biasIns_on = 0;

    % ### IMU OUTPUT ###
    IMUoutput = 'standard';

    % arm
    IMUarm = param.AC.IMUarm; %[0; 0; 0;]; %paramAC.IMUarm;
    wn = 20;

    param.sensors.IMU = struct('T_IMU',T_IMU,'IMUbias',IMUbias,'IMUarm',IMUarm,...
        'IMUmodel',IMUmodel,'IMUdrift',IMUdrift,'IMUtau',IMUtau,'wn',wn,...
        'Accstd',Accstd,'RateGyrostd',RateGyrostd,'IMUoutput',IMUoutput,...
        'IMU_DATAtype',IMU_DATAtype,'biasIns_on',biasIns_on);
end
T_IMU_print = param.sensors.IMU.T_IMU*1e3;
fprintf('IMU sample time: %3.2f ms \n',T_IMU_print)
fprintf('IMU output: %s \n',param.sensors.IMU.IMUoutput)

param.sim.IMUmodel = SetIMUmodel;
if SetIMUmodel ~= param.sensors.IMU.IMUmodel
    fprintf('IMU IMUmodel: %d => overided: %d \n',param.sensors.IMU.IMUmodel,...
        param.sim.IMUmodel)
else
    fprintf('IMU IMUmodel: %d \n',param.sensors.IMU.IMUmodel)
end

param.sim.IMUtau = SetIMUtau;
if SetIMUtau ~= param.sensors.IMU.IMUtau
    fprintf('IMU time constant: %d => overided: %d \n',param.sensors.IMU.IMUtau,...
        param.sim.IMUtau)
else
    fprintf('IMU time constant: %d s\n',param.sensors.IMU.IMUtau)
end

param.sim.biasIns_on = SetBiasIns_on;
if SetBiasIns_on ~= param.sensors.IMU.biasIns_on
    fprintf('IMU bias instability: %d => overided: %d \n',param.sensors.IMU.biasIns_on,...
        param.sim.biasIns_on)
else
    fprintf('IMU bias instability : %d \n',param.sensors.IMU.biasIns_on)
end

if strcmp(IMU_DATAtype,'SIM')
    fprintf('IMU type: %s \n\n',param.sensors.IMU.IMUtype)
elseif strcmp(IMU_DATAtype,'FT')
    fprintf('IMU type: %s \n\n',param.AC.IMUtype)
end

% ### PROCESS COVARIANCE ###
% FilterConfigParams
T_IMU       = param.sensors.IMU.T_IMU;
IMUbias     = param.sensors.IMU.IMUbias;
IMUdrift    = param.sensors.IMU.IMUdrift;
Accstd      = param.sensors.IMU.Accstd;
RateGyrostd = param.sensors.IMU.RateGyrostd;
% IMUtau      = param.sensors.IMU.IMUtau;
% Q matrix
stddevPos2      = [1 1 1];
stddevBias2     = Accstd^2/T_IMU*[1 1 1];
stddevDrift2    = RateGyrostd^2/T_IMU*[1 1 1];
if SetIMUmodel==2
    stddevBiasIns2  = [1 1 1]*(1e-3*g0)^2;
    stddevDriftIns2 = [1 1 1]*degph2radps^2;
elseif SetIMUmodel==3
    stddevBiasIns2  = 2*IMUbias^2/SetIMUtau*[1 1 1];
    stddevDriftIns2 = 2*IMUdrift^2/SetIMUtau*[1 1 1];
else
    stddevBiasIns2  = 0*[1 1 1];
    stddevDriftIns2 = 0*[1 1 1];
end
QIMUgain = diag([stddevPosGain.^2 stddevNablaGain.^2 stddevEpsilonGain.^2 ...
    stddevBiasInsGain.^2 stddevDriftInsGain.^2]);
Q_IMU = diag([stddevPos2 stddevBias2 stddevDrift2 stddevBiasIns2 stddevDriftIns2])*QIMUgain;
% QIMUgain = diag([stddevNablaGain.^2 stddevEpsilonGain.^2 ...
%     stddevBiasInsGain.^2 stddevDriftInsGain.^2]);
% Q_IMU = diag([stddevBias2 stddevDrift2 stddevBiasIns2 stddevDriftIns2])*QIMUgain;
% Initial continuous state covariance matrix
VarDR0 = [1 1 1];
VarDV0 = [1 1 1];
VarPsi0 = [1 1 1]*deg2rad^2;
VarNabla0 = IMUbias^2*[1 1 1];
VarEpsilon0 = IMUdrift^2*[1 1 1];
P0IMUgain = diag([stddevDR0Gain.^2 stddevDV0Gain.^2 stddevPsi0Gain.^2 ...
    stddevNabla0Gain.^2 stddevEpsilon0Gain.^2]);
P0_IMU = diag([VarDR0 VarDV0 VarPsi0 VarNabla0 VarEpsilon0])*P0IMUgain;
% ------------

% ### GT ###
disp('Loading Ground Truth data ...')
eval(['load GT_data_' FP '.mat'])
dt = tspan(2)-tspan(1);
dt = round(dt,6);
dt_print = dt*1e3;
fprintf('Ground Truth sample time: %4.3f ms \n\n',dt_print)
lat_GT  = X_GT(:,1)';
long_GT = X_GT(:,2)';
alt_GT  = X_GT(:,3)';
V_l_GT  = X_GT(:,4:6)';
q_bl_GT = X_GT(:,7:10)';

% ### SIMULATION SAMPLE TIME ###
switch param.sensors.IMU.IMUoutput
    case 'standard'
        T = round(param.sensors.IMU.T_IMU,6);
    case 'increments'
        T = round(4*param.sensors.IMU.T_IMU,6);
end
time = 0:T:tspan(end);
N = length(time);
ratioTdt = round(T/dt);
ratioTTIMU = round(T/T_IMU);
Tprint = T*1e3;
fprintf('Simulation sample time: %3.1f ms \n\n',Tprint)
if mod(T,dt)~=0
    warning('Simulation sample time is not multiple of Ground Truth sample time')
end
% param.sim = struct('T',T,'dt',dt,'MCreal',M);
param.sim.T = T;
param.sim.dt = dt;
param.sim.MCreal = M;
    
% ### PVA ###
if PVA_enabled
    disp('Loading PVA data ...')
    eval(['load Sensors/PVA/PVA_data_' FP '.mat'])
    PVA_DATAtype = 'SIM';
    param.sensors.PVA = paramPVA;
    clear paramPVA
    param.sensors.PVA.PVAmodel = SetPVAmodel;
    kPVA  = fix(param.sensors.PVA.T_PVA/T);
    
    T_PVA_print = param.sensors.PVA.T_PVA*1e3;
    fprintf('ALT sample time: %3.2f ms \n',T_PVA_print)
    % Tunning gain
    param.sensors.ALT.stdPVAposgain = stdPVAposgain;
    fprintf('std PVA Pos gain: %2.1f \n\n',stdPVAposgain)
    param.sensors.ALT.stdPVAvelgain = stdPVAvelgain;
    fprintf('std PVA Vel gain: %2.1f \n\n',stdPVAvelgain)
end

% ### GPS ###
if GPSTC_enabled
    % Load data
    disp('Loading GPS TC data ...')
    eval(['load Sensors/GPS/GPSTC_data_' FP '.mat'])
    % Load param
    if strcmp(GPSTC_DATAtype,'SIM')
        param.sensors.GPS = paramGPS;
        param.AC.GPSarm = paramGPS.GPSarm;
        clear paramGPS
    elseif strcmp(GPSTC_DATAtype,'FT')
        T_GPSTC = GPSTC_time(2)-GPSTC_time(1);
        
        GPSphasePSD  = 1.1e-19*c^2; % [m^2/s^2/Hz] - GPS clock phase PSD
        GPSfreqPSD   = 4.3e-20*c^2; % [m^2/s^4/Hz] - GPS clock frequency PSD
        GPSnoise_rho = 20;%5;%10/3;  % [m^2/s^4/Hz] - GPS rho white noise
        GPSnoise_vel = 2;%1;           % [m^2/s^4/Hz] - GPS rho white noise
        
        GPSarm = param.AC.GPSarm; %[0; 0; 0;]; %paramAC.GPSarm;
        
        param.sensors.GPS = struct('T_GPSTC',T_GPSTC,'GPSphasePSD',GPSphasePSD,...
            'GPSfreqPSD',GPSfreqPSD,'GPSnoise_rho',GPSnoise_rho,...
            'GPSnoise_vel',GPSnoise_vel,'GPSarm',GPSarm,'GPSTC_DATAtype',GPSTC_DATAtype,...
            'GPSFailureType',GPSFailureType,'Finterval',Finterval);
    end
    kGPSTC = fix(param.sensors.GPS.T_GPSTC/T);
    
    % GPS Failure Type
    if strcmp(param.sensors.GPS.GPSFailureType,'LOSS')
        Finterval = round(param.sensors.GPS.Finterval/60,1);
        fprintf('LOSS data between %d min and %d min \n',Finterval(1),Finterval(2))
    elseif strcmp(param.sensors.GPS.GPSFailureType,'ERRONEOUS')
        Finterval = round(param.sensors.GPS.Finterval/60,1);
        fprintf('ERRONEOUS data between %d min and %d min \n',Finterval(1),Finterval(2))
    end
    
    % Bias state
    nGPS = 2; % Number of GPS states
    P0_GPS = diag([stdClockBias0^2 stdClockBiasDot0^2]);
    QGPSgain = diag([QClockBiasGain^2 QClockBiasDotGain^2]);
    Q_GPS = diag([param.sensors.GPS.GPSphasePSD param.sensors.GPS.GPSfreqPSD])*QGPSgain; %<====
    
    % Ionospheric correction parameters
    alpha = [0.9313e-8 0 -0.5960e-7 0];
    beta  = [0.9011e5  0 -0.1966e6  0];
    
    param.sensors.GPS.ion_corr = struct('alpha',alpha,'beta',beta);
    
    % Trophosphere correction parameters
    P0      = [1013.25 1017.25 1015.75 1011.75 1013.00]';
    T0      = [299.65 294.15 283.15 272.15 263.65]';
    e0      = [26.31 21.79 11.66 6.78 4.11]';
    beta0   = [6.30e-3 6.05e-3 5.58e-3 5.39e-3 4.53e-3]';
    lambda0 = [2.77 3.15 2.57 1.81 1.55]';
    % Table0 = table(P0,T0,e0,beta0,lambda0);
    Table0 = [P0 T0 e0 beta0 lambda0];
    
    DP      = [0.00 -3.75 -2.25 -1.75 -0.50]';
    DT      = [0.00 7.00 11.00 15.00 14.50]';
    De      = [0.00 8.85 7.24 5.36 3.39]';
    Dbeta   = [0.00e-3 0.25e-3 0.32e-3 0.81e-3 0.62e3]';
    Dlambda = [0.00 0.33 0.46 0.74 0.30]';
    % TableDelta = table(DP,DT,De,Dbeta,Dlambda);
    TableDelta = [DP DT De Dbeta Dlambda];
    
    param.sensors.GPS.tropo_corr = struct('Table0',Table0,'TableDelta',TableDelta);
    
    % Receiver noise
    d = 0.2; % chip (correlator spacing)
    B = 0.5; % Hz (one sided code loop bandwidth)
    
    param.sensors.GPS.noise = struct('d',d,'B',B);
    
    T_GPSTC_print = param.sensors.GPS.T_GPSTC*1e3;
    fprintf('GPS TC sample time: %3.2f ms \n',T_GPSTC_print)
    % PR Tunning gain
    param.sensors.GPS.stdPRgain = stdPRgain;
    fprintf('std PR gain: %2.1f \n',stdPRgain)
    % DPR on
    param.sensors.GPS.DPRon = DPRon;
    fprintf('Delta range: %d \n',DPRon)
    % DPR Tunning gain
    param.sensors.GPS.stdDPRgain = stdDPRgain;
    fprintf('std DPR gain: %2.1f \n',stdDPRgain)
    % Earth rotation correction
    param.sensors.GPS.Erot_corr_enabled = Erot_corr_enabled;
    fprintf('Earth rotation correction: %d \n',Erot_corr_enabled)
    % Atm correction
    param.sensors.GPS.atm_corr_enabled = atm_corr_enabled;
    fprintf('Tropo/iono corrections: %d \n\n',atm_corr_enabled)
end
param.sensors.GPS.GPSTC_enabled = GPSTC_enabled;


% ### GPSLC ###
if GPSLC_enabled
    % Load data
    disp('Loading GPS LC data ...')
    eval(['load Sensors/GPS/GPSLC_data_' FP '.mat'])
    % Load param
    GPSLC_DATAtype = 'FT';
    if strcmp(GPSLC_DATAtype,'SIM')
        param.sensors.GPSLC = paramGPSLC;
        clear paramGPSLC
    elseif strcmp(GPSLC_DATAtype,'FT')
        param.sensors.GPSLC.T_GPSLC = GPSLC_time(2)-GPSLC_time(1);
        param.sensors.GPSLC.stdGPSLCpos = 5; 
        param.sensors.GPSLC.stdGPSLCvel = 1;
    end
    kGPSLC  = fix(param.sensors.GPSLC.T_GPSLC/T);
    T_GPSLC_print = param.sensors.GPSLC.T_GPSLC*1e3;
    fprintf('GPSLC sample time: %3.2f ms \n',T_GPSLC_print)
    % Tunning gain
    param.sensors.GPSLC.stdGPSLCgain = stdGPSLCgain;
    fprintf('std GPSLC gain: %2.1f \n\n',stdGPSLCgain)
end

% ### ALT ###
if ALT_enabled
    % Load data
    disp('Loading Altimeter data ...')
    eval(['load Sensors/ALT/ALT_data_' FP '.mat'])
    % Load param
    if strcmp(ALT_DATAtype,'SIM')
        param.sensors.ALT = paramALT;
        clear paramALT
    elseif strcmp(ALT_DATAtype,'FT')
        param.sensors.ALT.T_ALT = ALT_time(2) - ALT_time(1);
        param.sensors.ALT.stdALT = 10*ft2m;
    end
    kALT  = fix(param.sensors.ALT.T_ALT/T);
    % Bias state
    nALT = 1; % Number of ALT states
    Q_ALT = stdALTbias^2; % constant bias
    P0_ALT = stdALTP0^2; 

    T_ALT_print = param.sensors.ALT.T_ALT*1e3;
    fprintf('ALT sample time: %3.2f ms \n',T_ALT_print)
    % Tunning gain
    param.sensors.ALT.stdALTgain = stdALTgain;
    fprintf('std ALT gain: %2.1f \n\n',stdALTgain)
end
param.sensors.ALT.ALT_enabled = ALT_enabled;

% ### HDG ###
if HDG_enabled
    % Load data
    disp('Loading Heading data ...')
    eval(['load Sensors/HDG/HDG_data_' FP '.mat'])
    % Load param
    % HDG_DATAtype = 'FT';
    if strcmp(HDG_DATAtype,'SIM')
        param.sensors.HDG = paramHDG;
        clear paramHDG
    elseif strcmp(HDG_DATAtype,'FT')
        param.sensors.HDG.T_HDG = HDG_time(2)-HDG_time(1);
        param.sensors.HDG.stdHDG = 0.5*deg2rad;%2*deg2rad;
    end
    kHDG  = fix(param.sensors.HDG.T_HDG/T);
    T_HDG_print = param.sensors.HDG.T_HDG*1e3;
    fprintf('HDG sample time: %3.2f ms \n',T_HDG_print)
    % Tunning gain
    param.sensors.HDG.stdHDGgain = stdHDGgain;
    fprintf('std HDG gain: %2.1f \n\n',stdHDGgain)
end

% ### MAG ###
if MAG_enabled
    % Load data
    disp('Loading Magnetometer data ...')
    eval(['load Sensors/MAG/MAG_data_' FP '.mat'])
    % Load param
    MAG_DATAtype = 'SIM';
    param.sensors.MAG = paramMAG;
    clear paramMAG
    kMAG  = fix(param.sensors.MAG.T_MAG/T);
    
    T_MAG_print = param.sensors.MAG.T_MAG*1e3;
    fprintf('MAG sample time: %3.2f ms \n',T_MAG_print)
    % Tunning gain
    param.sensors.MAG.stdMAGgain = stdMAGgain;
    fprintf('std MAG gain: %2.1f \n\n',stdMAGgain)
end

% ### CAMLMF ###
if CAMLMF_enabled
    % Load data
    disp('Loading Pointing Camera Features Tracking data ...')
    eval(['load Sensors/CAM/CAMLMF/CAMLMF_data_' FP '.mat'])
    % Load param
    param.sensors.CAMLMF = paramCAMLMF;
    clear paramCAMLMF
    kCAMLMF = fix(param.sensors.CAMLMF.T_CAMLMF/T);
    T_CAMLMF_print = param.sensors.CAMLMF.T_CAMLMF*1e3;
    fprintf('CAM LMF sample time: %3.2f ms \n\n',T_CAMLMF_print)
    nLM = size(LM_time,1);
    fprintf('%d Landmarkds on the path. Visualization at following time intervals\n',nLM)
    for i=1:nLM
        fprintf('[%d %d] \n',LM_time(i,1),LM_time(i,2))
    end
    % Tunning gain
    param.sensors.CAMLMF.stdCAMgain = stdCAMgain;
    fprintf('stdCAMgain: %2.1f \n\n',stdCAMgain)
end

% ### CAMLMS ###
if CAMLMS_enabled
    % Load data
    disp('Loading Pointing Camera with Scale data ...')
    eval(['load Sensors/CAM/CAMLMS/CAMLMS_data_' FP '.mat'])
    % Load param
    % CAMLMF_DATAtype = 'SIM';
    param.sensors.CAMLMS = paramCAMLMS;
    param.sensors.CAMLMDB = paramCAMLMDB;
    clear paramCAMLMS paramCAMLMDB
    kCAMLMS = fix(param.sensors.CAMLMS.T_CAMLMS/T);
    T_CAMLMS_print = param.sensors.CAMLMS.T_CAMLMS*1e3;
    fprintf('CAM LMS sample time: %3.2f ms \n\n',T_CAMLMS_print)
    nLM = size(LM_time,1);
    fprintf('%d Landmarkds on the path. Visualization at following time intervals\n',nLM)
    for i=1:nLM
        fprintf('[%d %d] \n',LM_time(i,1),LM_time(i,2))
    end
    % Tunning gain
    param.sensors.CAMLMS.stdCAMgain = stdCAMgain;
    fprintf('stdCAMgain: %2.1f \n\n',stdCAMgain)
end

% ### CAMLMP ###
if CAMLMP_enabled
    % Load data
    disp('Loading Pointing Camera data ...')
    eval(['load Sensors/CAM/CAMLMP/CAMLMP_data_' FP '.mat'])
    % Load param
    CAMLMP_DATAtype = 'SIM';
    if strcmp(CAMLMP_DATAtype,'SIM')
        param.sensors.CAMLMP = paramCAMLMP;
        clear paramCAMLMP  
    elseif strcmp(CAMLMP_DATAtype,'FT')
        param.sensors.CAMLMP.T_CAMLMP = CAMLMP_time(2)-CAMLMP_time(1);
        param.sensors.CAMLMP.stdCAMLMP = 2*deg2rad;
    end
    kCAMLMP = fix(param.sensors.CAMLMP.T_CAMLMP/T);
    T_CAMLMP_print = param.sensors.CAMLMP.T_CAMLMP*1e3;
    fprintf('CAM LMP sample time: %3.2f ms \n\n',T_CAMLMP_print)
    nLM = size(LM_time,1);
    fprintf('%d Landmarkds on the path. Visualization at following time intervals\n',nLM)
    for i=1:nLM
        fprintf('[%d %d] \n',LM_time(i,1),LM_time(i,2))
    end
    % Tunning gain
    param.sensors.CAMLMP.stdCAMgain = stdCAMgain;
    fprintf('stdCAMgain: %2.1f \n\n',stdCAMgain)
end

% ### CAMLM ###
if CAMLM_enabled
    % Load data
    disp('Loading Stabilized Camera data ...')
    eval(['load Sensors/CAM/CAMLM/CAMLM_data_' FP '.mat'])
    % Load param
    CAMLM_DATAtype = 'SIM';
    if strcmp(CAMLM_DATAtype,'SIM')
        param.CAMLM = paramCAMLM;
        clear paramCAMLM
    elseif strcmp(CAMLM_DATAtype,'FT')
        param.CAMLM.T_CAMLM = CAMLM_time(2)-CAMLM_time(1);
        param.CAMLM.stdCAMLM = 2*deg2rad;
    end
    kCAMLM = fix(param.CAMLM.T_CAMLM/T);
end

% ### LPS ###
if LPS_enabled
    % Load data
    disp('Loading LPS data ...')
    eval(['load Sensors/LPS/LPS_data_' FP '.mat'])
    % Load param
    param.sensors.LPS = paramLPS;
    clear paramLPS
    kLPS = fix(param.sensors.LPS.T_LPS/T);

    T_LPS_print = param.sensors.LPS.T_LPS*1e3;
    fprintf('LPS sample time: %3.2f ms \n',T_LPS_print)
    % Tunning gain
    param.sensors.LPS.stdLPSgain = stdLPSgain;
    fprintf('std LPS gain: %2.1f \n\n',stdLPSgain)
end

% -------------------------------------------------------------------------

% ### CAMOF ###
if CAMOF_enabled
    % Load data
    disp('Loading Opticial Flow Camera data ...')
    load Sensors/CAMOF/CAMOF_data.mat
    
    T_CAMOF = time_CAMOF(2)-time_CAMOF(1);
    kCAMOF = fix(T_CAMOF/T);
    
    fOF = 20/(32/512);%100;
    pan_OF = 0;
    tilt_OF = -90*deg2rad;
    stdCAMOF = 1.5*10;
    
    param.CAMOF = struct('f',fOF,'pan_OF',pan_OF,'tilt_OF',tilt_OF,...
        'stdCAMOF',stdCAMOF);
end

% ### SR ###
if SR_enabled
    % Load data
    disp('Loading Slant Range data ...')
    eval(['load Sensors/SlantRange/SR_data_' FP '.mat'])
    % Load param
    param.SR = paramSR;
    clear paramSR
    kSR = fix(param.SR.T_SR/T);
end

% ### WOW ###
WOW_enabled = false;