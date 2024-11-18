close all
clear 
clc 

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

%% Script settings
%--------------------------------------------------------------------------
subidx  = 1;
subject = ['sub-',num2str(subidx,'%02d')];

datapath = cell(1,2);
datapath{1} = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-upchirps.fif']);
datapath{2} = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-upchirps_fail.fif']);

avg = cell(1,2);

%% Preprocess data
%--------------------------------------------------------------------------

for i=1:2
    % Define trials
    %--------------
    cfg                     = [];
    cfg.dataset             = datapath{i};
    cfg.trialfun            = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = 2; 
    cfg.trialdef.prestim    = 0.05;                  % in seconds
    cfg.trialdef.poststim   = 0.35;                  % in seconds
    cfg                     = ft_definetrial(cfg);
    trl                     = cfg.trl;
                
    % Filter continuous data to avoid edge artifacts
    %-----------------------------------------------
    cfg              = [];
    cfg.dataset      = datapath{i};
    cfg.channel      = 'meg'; 
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
    data    = ft_redefinetrial(cfg,data); 

    % Demean trials
    %--------------
    cfg                = [];
    cfg.demean         = 'yes';
    cfg.baselinewindow = [-0.05 0];
    data               = ft_preprocessing(cfg,data);  
    
    %% Timelockanalysis
    cfg    = [];
    avg{i} = ft_timelockanalysis(cfg, data);

end % loop conditions

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
% cfg.layout     = 'neuromag306planar.lay';
cfg.linecolor  = 'rb';
ft_multiplotER(cfg, avg{:});
legend({'upchirps';'upchirps_fail'});



%% browse data

% Define trials
%--------------
cfg                     = [];
cfg.dataset             = datapath{1};
cfg.trialfun            = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype  = 'STI101';
cfg.trialdef.eventvalue = 2; 
cfg.trialdef.prestim    = 0.05;                  % in seconds
cfg.trialdef.poststim   = 0.35;                  % in seconds
cfg                     = ft_definetrial(cfg);
trl                     = cfg.trl;

% Filter continuous data to avoid edge artifacts
%-----------------------------------------------
cfg              = [];
cfg.dataset      = datapath{i};
cfg.channel      = 'meg'; 
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
data    = ft_redefinetrial(cfg,data); 

cfg          = [];
cfg.channel  = {'MEGMAG','-MEG1813','-MEG2622'};
cfg.method   = 'summary';
dummy        = ft_rejectvisual(cfg, data);

cfg          = [];
cfg.channel  = 'MEGGRAD';
cfg.method   = 'summary';
dummy        = ft_rejectvisual(cfg, data);
    
%% Timelockanalysis
cfg       = [];
avg_dummy = ft_timelockanalysis(cfg, dummy);

cfg      = [];
avg_data = ft_timelockanalysis(cfg, data);

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
% cfg.layout     = 'neuromag306planar.lay';
cfg.linecolor  = 'rb';
ft_multiplotER(cfg,avg_data,avg_dummy);
legend({'data';'data cleaned'});

%% databrowser

% Define trials
%--------------
cfg                     = [];
cfg.dataset             = datapath{1};
cfg.trialfun            = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype  = 'STI101';
cfg.trialdef.eventvalue = 2; 
cfg.trialdef.prestim    = 0.05;                  % in seconds
cfg.trialdef.poststim   = 0.35;                  % in seconds
cfg                     = ft_definetrial(cfg);
trl                     = cfg.trl;

% Filter continuous data to avoid edge artifacts
%-----------------------------------------------
cfg              = [];
cfg.dataset      = datapath{1};
cfg.channel      = 'meg'; 
cfg.continuous   = 'yes';
cfg.bpfilter     = 'yes';
cfg.bpfreq       = [16,120];
cfg.dftfilter    = 'yes';       % enable notch filtering to eliminate power line noise
cfg.dftfreq      = [50,100];    % set up the frequencies for notch filtering
cfg.dftreplace   = 'zero';      % 'zero' is default'
cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
data             = ft_preprocessing(cfg); 

% ecg artifacts
cfg                 = [];
cfg.continuous      = 'yes';
cfg.dataset         = datapath{1};
cfg.artfctdef.ecg.channel = {'ECG003'};
cfg.channel         = {'ECG003'};   
[~, artifact_ecg]   = ft_artifact_ecg(cfg);

cfg           = [];
cfg.blocksize = 10;
cfg.viewmode  = 'butterfly';
cfg.channel   = 'megmag';
% cfg.trl       = trl;
cfg.artfctdef.xxx.artifact
artf          = ft_databrowser(cfg, data);

%% ICA test











rmpath('..')


