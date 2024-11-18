close all
clear 
clc 

% Setup
% sub-01: 42 channel EEG, 32 channel EEG cap + 10 extra electrodes around ear
% sub-02: 10 channel ear EEG

%% Script settings
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

subidx  = 2; 
subject = ['sub-',num2str(subidx,'%02d')];

conditions = {'clicks','upchirps','downchirps'};

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

% Directory for data storing
%---------------------------
dir2save = fullfile(settings.path2project,'derivatives',subject,'sensorlevel');

% Channel names 
%--------------
switch subject
    case 'sub-01'
        ear_channels = {'EEG038','EEG039','EEG040','EEG041','EEG042','EEG008',... % left
                        'EEG033','EEG034','EEG035','EEG036','EEG037','EEG012'}; % right
    case 'sub-02'
        ear_channels = {'EEG006','EEG007','EEG008','EEG009','EEG010',... % left
                        'EEG001','EEG002','EEG003','EEG004','EEG005'}; % right
end

%% Preprocess data - Epoching
%--------------------------------------------------------------------------
epochs_condition = cell(1,3);
N_trials         = zeros(1,3);

for cidx = 1:3 % loop over conditions

    % Path for EEG data
    %------------------
    condition = conditions{cidx};
    datapath  = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);

    % Browse data
    %------------
    % cfg           = [];
    % cfg.dataset   = datapath;
    % cfg.channel   = ear_channels;
    % cfg.blocksize = 10;
    % cfg.preproc.detrend = 'yes';
    % ft_databrowser(cfg);

    switch condition
        case 'clicks'
            eventvalue = TrigID.click;
        case 'upchirps'
            eventvalue = TrigID.upchirp;
        case 'downchirps'
            eventvalue = TrigID.downchirp;
    end

    % Define trials
    %--------------
    cfg                     = [];
    cfg.dataset             = datapath;
    cfg.trialfun            = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = eventvalue; 
    cfg.trialdef.prestim    = 0.05;                  % in seconds
    cfg.trialdef.poststim   = 0.35;                  % in seconds
    cfg                     = ft_definetrial(cfg);
    trl                     = cfg.trl;

    if ~isequal(length(trl),1200)
        error('Unexpected number of triggers.')
    end

    % Filter continuous data to avoid edge artifacts
    %-----------------------------------------------
    cfg              = [];
    cfg.dataset      = datapath;
    cfg.channel      = ear_channels;
    cfg.continuous   = 'yes';
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = [16,120];
    cfg.dftfilter    = 'yes';       % enable notch filtering to eliminate power line noise
    cfg.dftfreq      = [50,100];    % set up the frequencies for notch filtering
    cfg.dftreplace   = 'zero';      % 'zero' is default'
    cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
    data             = ft_preprocessing(cfg); 

    % Epoch data
    %-----------
    cfg     = [];
    cfg.trl = trl;  
    epochs  = ft_redefinetrial(cfg,data); 

    % Demean trials
    %--------------
    % cfg                = [];
    % cfg.demean         = 'yes';
    % cfg.baselinewindow = [-0.05 0];
    % epochs            = ft_preprocessing(cfg,epochs);  

    % Semi automatic artifact rejection
    %----------------------------------
    cfg          = [];
    cfg.metric   = 'maxzvalue';
    epochs_clean = ft_rejectvisual(cfg,epochs);

    N_trials(cidx) = length(epochs_clean.trial);
    epochs_condition{cidx} = epochs_clean;

    clear data epochs epochs_clean
end % conditions

%% Save epochs
%--------------------------------------------------------------------------
results.epochs     = epochs_condition; 
results.conditions = conditions;
results.N_trials   = N_trials;

save(fullfile(dir2save,[subject,'_chirps_erps_earEEG.mat']),'-struct','results'); 