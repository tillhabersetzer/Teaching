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

conditions = settings.conditions;

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

% ICA components
components         = importdata(fullfile(settings.path2project,'derivatives',subject,'sensorlevel',[subject,'_ica.mat']));
bad_components_meg = {'sub-01',[2,21];...
                      'sub-02',[]};

data_preprocessed     = cell(1,3);
avg                   = cell(1,3);
data_preprocessed_ica = cell(1,3);
avg_ica               = cell(1,3);

for cidx=1:3 % loop over conditions

    % Path for MEG data
    %--------------------------------------------------------------
    % Path for MEG data (no headposition transformation applied)
    if settings.maxfilter
        datapath = fullfile(settings.path2project,'derivatives',subject,'meg',[subject,'_task-',conditions{cidx},'_tsss.fif']);
    else
        datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',conditions{cidx},'.fif']);
    end

    switch conditions{cidx}
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
            
    % Filter continuous data to avoid edge artifacts
    %-----------------------------------------------
    cfg              = [];
    cfg.dataset      = datapath;
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

    % Apply ICA
    %----------
    cfg           = [];
    idx           = find(contains(bad_components_meg(:,1),subject));
    cfg.component = bad_components_meg{idx,2};
    data_ica      = ft_rejectcomponent(cfg, components, data);

    % Demean trials
    %--------------
    cfg                = [];
    cfg.demean         = 'yes';
    cfg.baselinewindow = [-0.05 0];
    data               = ft_preprocessing(cfg,data);  
    data_ica           = ft_preprocessing(cfg,data_ica);  

    % Rename data for pooling
    %------------------------
    data_preprocessed{cidx}     = data;
    data_preprocessed_ica{cidx} = data_ica;

    %% Timelockanalysis
    cfg           = [];
    avg{cidx}     = ft_timelockanalysis(cfg, data);
    avg_ica{cidx} = ft_timelockanalysis(cfg, data_ica);

end % loop conditions

%% comparison

layout = 'neuromag306mag.lay';
% layout = 'neuromag306planar.lay';

% all conditions
%---------------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbk';
cfg.showlegend = 'yes' ;
ft_multiplotER(cfg, avg{:});
legend({'clicks';'upchirps';'downchirps'});
subtitle(['without ICA | ' subject])

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbk';
ft_multiplotER(cfg, avg_ica{:});
legend({'clicks';'upchirps';'downchirps'});
title(['with ICA | ' subject])

% compare conditions
%-------------------
cidx = 2;

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
ft_multiplotER(cfg, avg{cidx},avg_ica{cidx});
legend({'without ICA','with ICA'});
title([conditions{cidx} '|' subject])



cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
ft_multiplotER(cfg, avg{cidx});
title([conditions{cidx} '|' subject])



