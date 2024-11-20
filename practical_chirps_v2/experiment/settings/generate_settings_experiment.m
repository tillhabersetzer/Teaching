%--------------------------------------------------------------------------
% Till Habersetzer, 15.10.2024
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Revision:
% 05.11.24
% 15.11.24
%--------------------------------------------------------------------------

close all
clearvars
clc

%% Calibration: Specify values in dB (p-p) peSPL
%--------------------------------------------------------------------------

% Calibration Values for different earphones
%-------------------------------------------

% TIP 300 earphones without equalization
settings.calibration.tip300.cal_val.click = [111.0,111.5];
settings.calibration.tip300.cal_val.up    = [116.5,116.5];
settings.calibration.tip300.cal_val.down  = [116.5,116.5];

% Sensimetrics earphones without equalization
settings.calibration.sensimetrics.cal_val.click = [107.5,109.5];
settings.calibration.sensimetrics.cal_val.up    = [114.5,115.5];
settings.calibration.sensimetrics.cal_val.down  = [114,116];

% Sensimetrics earphones with equalization
settings.calibration.sensimetrics_eq.cal_val.click = [102,104];
settings.calibration.sensimetrics_eq.cal_val.up    = [106.5,110];
settings.calibration.sensimetrics_eq.cal_val.down  = [106.5,111];

%% Threshod estimates in dB (p-p) peSPL 
%--------------------------------------------------------------------------

% Start Level in adaptive hearing procedure
%------------------------------------------

% Dau 2000
% mean(0 dB SL, click) ~ 45 dB peSPL
% mean(0 dB SL, up) ~ 33 dB peSPL

% Start 20 dB over dummy threshold for stimuli
settings.threshold.start_level.click  = 60;
settings.threshold.start_level.up     = 60;
settings.threshold.start_level.down   = 60;

% TIP 300 earphones without equalization
settings.threshold.tip300.click = [33,33]; 
settings.threshold.tip300.up    = [30,30]; 

% Sensimetrics earphones without equalization
settings.threshold.sensimetrics.click = [32,32];
settings.threshold.sensimetrics.up    = [32,32];

% Sensimetrics earphones with equalization
settings.threshold.sensimetrics_eq.click = [33,33]; 
settings.threshold.sensimetrics_eq.up    = [33,33]; 

% Other settings for calibration signal
settings.level_max = 100; % maximum level 100 dB (p-p) peSPL 

%% SoundMexPro Settings
%--------------------------------------------------------------------------
% settings.asio_driver = 'ASIO Fireface USB'; 
settings.asio_driver = 'Focusrite USB ASIO';
settings.samp_freq   = 44100;

%% Settings for experiment
%--------------------------------------------------------------------------

% Define Trigger IDs
%-------------------
settings.trig_id.click = 1;
settings.trig_id.up    = 2;
settings.trig_id.down  = 4;

% Others
%--------------------------------------------------------------------------
settings.setup              = 'Practical MEG';
settings.n_stimuli          = 1200;
settings.jitter             = [0.35,0.05];  % 350 ms + 50 ms randomized jitter
settings.target_level_dB_SL = 40;   

% Save Settings
%--------------
settings.last_update = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % time of last update            
save('settings_experiment.mat','settings')