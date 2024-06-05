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

% only valid for subject 01: contains cap
subidx  = 1; 
subject = ['sub-',num2str(subidx,'%02d')];

conditions = {'clicks','upchirps','downchirps'};
%--------------------------------------------------------------------------

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

%% Create EEG layout
datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-clicks.fif']);
elec     = ft_read_sens(datapath, 'senstype', 'eeg');

ft_plot_sens(elec,'label','on','elecshape','disc');

cfg         = [];
cfg.elec    = elec;
cfg.channel = 'all'; 
layout      = ft_prepare_layout(cfg);

cfg        = [];
cfg.layout = layout;
ft_layoutplot(cfg)

%% Preprocess data 
%--------------------------------------------------------------------------
avgs     = cell(1,3);
N_trials = zeros(1,3);

for cidx = 1:3 % loop over conditions

    % Path for EEG data
    %------------------
    condition = conditions{cidx};
    datapath  = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);

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
    cfg.channel      = 'eeg';
    cfg.continuous   = 'yes';
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = [16,120];
    cfg.dftfilter    = 'yes';       % enable notch filtering to eliminate power line noise
    cfg.dftfreq      = [50,100];    % set up the frequencies for notch filtering
    cfg.dftreplace   = 'zero';      % 'zero' is default'
    cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
    data             = ft_preprocessing(cfg); 

    % Re-reference
    %-------------
    cfg            = [];
    cfg.reref      = 'yes';
    cfg.refmethod  = 'avg';
    cfg.channel    = 'all'; % this is the default
    cfg.refchannel = 'all'; % average is computed over the channels specified in cfg.refchannel
    data           = ft_preprocessing(cfg,data);

    % Epoch data
    %-----------
    cfg     = [];
    cfg.trl = trl;  
    epochs  = ft_redefinetrial(cfg,data); 

    % Demean trials
    %--------------
    cfg                = [];
    cfg.demean         = 'yes';
    cfg.baselinewindow = [-0.05 0];
    epochs            = ft_preprocessing(cfg,epochs);  

    % Semi automatic artifact rejection
    %----------------------------------
    cfg          = [];
    cfg.metric   = 'maxzvalue';
    epochs_clean = ft_rejectvisual(cfg,epochs);

    % Timelockanalysis
    %-----------------
    cfg        = [];
    avgs{cidx} = ft_timelockanalysis(cfg, epochs_clean);

    % Note down amount of preserved epochs
    %-------------------------------------
    N_trials(cidx) = length(avgs{cidx}.cfg.previous.trials);

    clear data epochs epochs_clean
end % conditions

%% Visualize data

% topoplot
%---------
cfg            = [];
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.showlabels = 'yes';
cfg.linecolor  = 'rbg';
ft_multiplotER(cfg, avgs{:});
legend(conditions);
subtitle('avg reference')

% single channel
%---------------
chan2plot = 'EEG005';
cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avgs{:});
legend(conditions);

%% Cleaning
rmpath('..')
