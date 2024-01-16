%% =========================== MNS SETTINGS ===============================

% ### SELECT FLIGHT PROFILE ###
% FP4
% FTE2
FP = 'FP4';
% ----------------------

% ### SENSOR LIST ###
sensorsList = {'IMU','WOW','PVA','ALT','GPSTC','GPSLC','MAG','HDG',...
    'CAMLMF','CAMLMS','CAMLMP','CAMLM','CAMOF','SR','LPS'};

% ### MAIN SENSOR => KALMAN FILTER PREDICTION ###
% IMU: Inertial Measurement Unit (accelerometers ans gyrometers) -> always
% enabled

% ### AUXILIARY SENSORS => KALMAN FILTER UPDATE ###
% WOW: Weight on wheels sensor
% PVA: Position, velocity and attitude observation
% ALT: Barometric altimeter
% GPSTC: GPS with pseudo range / tightly coupled
% GPSLC: GPS with LLA / loosely coupled
% MAG: Magnetometer / Earth magnetic field measurement
% HDG: Magnetic heading / compass
% CAMLMP: Landmark localization by tilt-pan pointing camera
% CAMLM: Landmark localization by tilt-pan stabilized camera
% CAMOF: Optical Flow by tilt-pan camera
% SR: Slant Range sensor
% LPS1: Landing Precision System Vertport 1
% LPS2: Landing Precision System Vertport 2

% ### SELECTED AUXILIARY SENSOR ###
sensorSelected = {'GPSTC','MAG','LPS'};%{'GPSTC','ALT','MAG','CAMLMF','LPS'};
% sensorSelected = {'PVA','CAMLMF'};
% ----------------------

% ### SELECT SENSOR DATA TYPE ###
% DEFAULT: SIM
% For FT data, select put the desired sensors in the list
FTList = {};
%FTList ={'IMU','ALT','GPSTC'};
% ----------------------

% ### KALMAN FILTER SETTINGS ###
% Filter Type
filter = 'EKF'; %'LKF' %'EKF

starWithFeedback = true;
tracePnorm_threshold = 25;%1e10;%1.0;
correction_time = 0;%8*60;%120;
TSP_start = 0;%8*60;%240;

% ChiSquare Test
chi2per = 0.95;
ChiSquareTestOff = true;
ChiSquareGain = 1;
Underweighting = false;

TRIAD_on = false;
TRIAD_MAG = true;

% MC realizations
M = 20; 

param.config = struct('FP',FP,'tracePnorm_threshold',tracePnorm_threshold,...
    'correction_time',correction_time,'TSP_start',TSP_start,...
    'chi2per',chi2per,'ChiSquareTestOff',ChiSquareTestOff,'MCreal',M,...
    'ChiSquareGain',ChiSquareGain,'Underweighting',Underweighting);
% -----------
saveData = true;

% -------------------------- SENSORS CONFIG -------------------------------
% IMU
% Load
IMUtype = 'HQTG';
IMUmodel = 1; %1-constant bias; 2- Bias Random Walk; 3- Gauss-Markov; 4- GM no noise;
% Set
SetIMUmodel = 2; %<===============
SetIMUtau = 3600;%1e12 %3600;%<===============
SetBiasIns_on = 1; %<===============
% Q matrix gain                         | FP4 - LQTG 1  ||  FP4 - LQTG MC ||    FTE2     |
stddevPosGain      = [1 1 1]*0.5;
stddevNablaGain    = 0.2*[1 1 2];           %| [1/5 1/5 1]/2 ||  [1 1 4]*5     ||  [1 1 1]     |
stddevEpsilonGain  = 0.2*[1 1 1];         %| [1 1 1/2]/5   ||  [1 1 2]*5    ||  [1 1 1]/3   |
stddevBiasInsGain  = 1e-4*[1 1 1];          %| [1 1 1]*10    ||  [1 1 1]      ||  [1 1 1]     |
stddevDriftInsGain = 1e-4*[1 1 1];          %| [1 1 1]*10    ||  [1 1 1]       |
% P0 matrix gain
stddevDR0Gain      = 20*[1 1 1];      %| 10*[1 1 1]     ||  25*[1 1 2]   ||  5*[1 1 1]   |
stddevDV0Gain      = [2 2 5];       %| 2*[1 1 1]      ||  5*[1 1 2]    ||  2*[1 1 1]   |
stddevPsi0Gain     = [2 2 4];       %| 1*[1 1 2]     ||  10*[1 1 1]    ||  2*[1 1 5]   |
stddevNabla0Gain   = 10*[1 1 1];       %| 2*[2 2 1]     ||  20*[1 1 1]    ||  [1 1 1]     |
stddevEpsilon0Gain = 10*[1 1 1];      %| 5*[1 1 1]     ||  20*[1 1 1]    ||  20*[1 1 1]  |
% ALT
stdALTbias = 0.01; % Bias process std   %|     0         ||     0         ||    0    |
stdALTP0 = 0.05;% Bias init cov         %|    0.25       ||    0.25       ||   0.15      |
stdALTgain = 4; % Tunning gain       %|     5         ||    1          ||    5    |
% GPSTC
stdClockBias0 = 10; % m              %|    5          ||    10        ||   10   |
stdClockBiasDot0 = 5;%50; % (m/s)   %|    100         ||    50        ||  150   |
QClockBiasGain = 1; %              %|    1/2        ||    1         ||  1/4   |
QClockBiasDotGain = 1;              %|    1          ||    1         ||  1/2   |
stdPRgain = 1; % Tunning gain       %|    5          ||    1         ||  1     |
stdDPRgain = 1; % Tunning gain      %|    5          ||    1         ||  1     |
DPRon = 1;                        %|    0          ||    0          ||  0     |  Delta range on/off
Erot_corr_enabled = 1;              %|    1          ||    1          ||  0     |  Earth rotation correction
atm_corr_enabled = 0;               %|    0          ||    0          ||  1     |  Atm correction
% GPSLC
stdGPSLCgain = 1;                   %|              ||    1          ||      |
% HDG
stdHDGgain = 1;                     %|    5          ||    1          ||   10   |
% MAG
stdMAGgain = 1;                     %|    5          ||    1          ||   10   |
% CAMLMP
% NOTES: stdCAMgain = 20 works for IMU FT + CAMLMP
stdCAMgain = 1;  % Tunning gain   %|    10         ||    2          ||   10   |
% LPS
stdLPSgain = 1;                    %|    5          ||    1          ||      |
% PVA
SetPVAmodel = 4;
stdPVAposgain = 1;
stdPVAvelgain = 1;

param.settings = struct('SetIMUmodel',SetIMUmodel,'SetIMUtau',SetIMUtau,...
    'SetBiasIns_on',SetBiasIns_on,...
    'stddevPosGain',stddevPosGain,'stddevNablaGain',stddevNablaGain,...
    'stddevEpsilonGain',stddevEpsilonGain,'stddevBiasInsGain',stddevBiasInsGain,...
    'stddevDriftInsGain',stddevDriftInsGain,'stddevDR0Gain',stddevDR0Gain,...
    'stddevDV0Gain',stddevDV0Gain,'stddevPsi0Gain',stddevPsi0Gain,...
    'stddevNabla0Gain',stddevNabla0Gain,'stddevEpsilon0Gain',stddevEpsilon0Gain,...
    'stdALTbias',stdALTbias,'stdALTP0',stdALTP0,'stdALTgain',stdALTgain,...
    'stdClockBias0',stdClockBias0,'stdClockBiasDot0',stdClockBiasDot0,...
    'QClockBiasGain',QClockBiasGain,'QClockBiasDotGain',QClockBiasDotGain,...
    'stdPRgain',stdPRgain,'stdDPRgain',stdDPRgain,'DPRon',DPRon,...
    'Erot_corr_enabled',Erot_corr_enabled,'atm_corr_enabled',atm_corr_enabled,...
    'stdGPSLCgain',stdGPSLCgain,'stdHDGgain',stdHDGgain,'stdMAGgain',stdMAGgain,...
    'stdCAMgain',stdCAMgain,'stdLPSgain',stdLPSgain,'SetPVAmodel',SetPVAmodel,...
    'stdPVAposgain',stdPVAposgain,'stdPVAvelgain',stdPVAvelgain);


