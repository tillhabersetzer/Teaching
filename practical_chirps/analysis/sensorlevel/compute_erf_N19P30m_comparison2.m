close all
clear 
clc 

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

% check whether maxfiltered data should be analyzed
maxfilter = settings.maxfilter;

%% Script settings
%--------------------------------------------------------------------------
subidx  = 2;
%--------------------------------------------------------------------------

% conditions
conditions = settings.conditions;

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

%% 1.) Preprocess data - Loose trial rejection
%--------------------------------------------------------------------------
subject           = ['sub-',num2str(subidx,'%02d')];
data_preprocessed = cell(1,3);
avg_loose         = cell(1,3);
N_trials_loose    = zeros(1,3);

for cidx=1:3 % loop over conditions

    % Path for MEG data
    %------------------------------------------------------------------
    condition = conditions{cidx};

    if maxfilter
        datapath = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-',condition,'-raw_tsss.fif']);
    else
        datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);
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
    cfg.trialfun            = 'ft_trialfun_general'; 
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
    epoched = ft_redefinetrial(cfg,data); 

    % Semi automatic artifact rejection
    %----------------------------------
    % separately for magnetometers
    cfg             = [];
    cfg.metric      = 'zvalue';
    cfg.channel     = 'megmag';
    cfg.keepchannel = 'yes';  % This keeps those channels that are not displayed in the data
    epoched         = ft_rejectvisual(cfg,epoched);
    
    % separately for gradiometers
    cfg.channel     = 'megplanar';
    epoched         = ft_rejectvisual(cfg,epoched);

    % Demean trials
    %--------------
    cfg                = [];
    cfg.demean         = 'yes';
    cfg.baselinewindow = [-0.05 0];
    epoched            = ft_preprocessing(cfg,epoched);  

    % Note down amount of preserved epochs
    %-------------------------------------
    N_trials_loose(cidx) = length(epoched.trial);

    % Timelockanalysis
    %-----------------
    cfg             = [];
    avg_loose{cidx} = ft_timelockanalysis(cfg, epoched);

    clear data epoched

end % loop over conditions

%% 2.) Preprocess data - Tight trial rejection
%--------------------------------------------------------------------------
subject           = ['sub-',num2str(subidx,'%02d')];
data_preprocessed = cell(1,3);
avg_tight         = cell(1,3);
N_trials_tight    = zeros(1,3);

for cidx=1:3 % loop over conditions

    % Path for MEG data
    %------------------------------------------------------------------
    condition = conditions{cidx};

    if maxfilter
        datapath = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-',condition,'-raw_tsss.fif']);
    else
        datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);
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
    cfg.trialfun            = 'ft_trialfun_general'; 
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
    epoched = ft_redefinetrial(cfg,data); 

    % Semi automatic artifact rejection
    %----------------------------------
    % separately for magnetometers
    cfg             = [];
    cfg.metric      = 'zvalue';
    cfg.channel     = 'megmag';
    cfg.keepchannel = 'yes';  % This keeps those channels that are not displayed in the data
    epoched         = ft_rejectvisual(cfg,epoched);
    
    % separately for gradiometers
    cfg.channel     = 'megplanar';
    epoched         = ft_rejectvisual(cfg,epoched);

    % Demean trials
    %--------------
    cfg                = [];
    cfg.demean         = 'yes';
    cfg.baselinewindow = [-0.05 0];
    epoched            = ft_preprocessing(cfg,epoched);  

    % Note down amount of preserved epochs
    %-------------------------------------
    N_trials_tight(cidx) = length(epoched.trial);

    % Timelockanalysis
    %-----------------
    cfg             = [];
    avg_tight{cidx} = ft_timelockanalysis(cfg, epoched);

    clear data epoched

end % loop over conditions

%% Compare results

%% Plot of event related fields for all sensors arranged topographically 
%--------------------------------------------------------------------------
layout = 'neuromag306mag.lay';
% layout = 'neuromag306planar.lay';

% without maxfilter
%------------------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbg';
ft_multiplotER(cfg, avg_loose{1}, avg_loose{2}, avg_loose{3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('loose')

% with maxfilter
%---------------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbg';
ft_multiplotER(cfg, avg_tight{1}, avg_tight{2}, avg_tight{3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('tight')

% comparison
%-----------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rb';
ft_multiplotER(cfg, avg_loose{1}, avg_tight{1});
legend({'loose','tight'})
subtitle(conditions{1});

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rb';
ft_multiplotER(cfg, avg_loose{2}, avg_tight{2});
legend({'loose','tight'})
subtitle(conditions{2});

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rb';
ft_multiplotER(cfg, avg_loose{3}, avg_tight{3});
legend({'loose','tight'})
subtitle(conditions{3});

%% Plot one specific channel
%--------------------------------------------------------------------------
% plot all 3 conditions together
chan2plot = 'MEG1411';

cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg_loose{1,:});
legend({conditions{1};conditions{2};conditions{3}});
title(['loose: ',subject])

cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg_tight{1,:});
legend({conditions{1};conditions{2};conditions{3}});
title(['tight: ',subject])

%% Visualize data - topography
timewin = [0.02 0.05];
layout = 'neuromag306mag.lay';
% layout = 'neuromag306planar.lay';

cfg                 = [];
cfg.layout          = layout;
cfg.xlim            = timewin;
cfg.style           = 'straight';
cfg.comment         = 'no';
cfg.marker          = 'off';
cfg.colorbar        = 'southoutside';

cfg.figure = subplot(2,3,1);
ft_topoplotER(cfg, avg_loose{1});
title(['loose: ',conditions{1}]);

cfg.figure = subplot(2,3,2);
ft_topoplotER(cfg, avg_loose{2});
title(['loose: ',conditions{2}]);

 cfg.figure= subplot(2,3,3);
ft_topoplotER(cfg, avg_loose{3});
title(['loose: ',conditions{3}]);

cfg.figure = subplot(2,3,4);
ft_topoplotER(cfg, avg_tight{3});
title(['tight: ',conditions{1}]);

cfg.figure = subplot(2,3,5);
ft_topoplotER(cfg, avg_tight{2});
title(['tight: ',conditions{2}]);

 cfg.figure= subplot(2,3,6);
ft_topoplotER(cfg, avg_tight{3});
title(['tight: ',conditions{3}]);

%% Combining planar gradiometers
% To help with identifying underlying sources we should make use of the 
% other channels in the data. The planar gradiometers are often more easily 
% interpreted, because they are most sensitive right above a source. 

% Combine planar
avg_loose_cmb = cell(1,3);
avg_tight_cmb = cell(1,3);
cfg     = [];
for n = 1:3
    avg_loose_cmb{n} = ft_combineplanar(cfg, avg_loose{n});
    avg_tight_cmb{n} = ft_combineplanar(cfg, avg_tight{n});
end

% Plot results
%-------------
cfg        = [];
cfg.layout = 'neuromag306cmb.lay';
ft_multiplotER(cfg, avg_loose_cmb{1}, avg_loose_cmb{2}, avg_loose_cmb{3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('loose')

cfg          = [];
cfg.layout   = 'neuromag306cmb.lay';
ft_multiplotER(cfg, avg_tight_cmb{1}, avg_tight_cmb{2}, avg_tight_cmb{3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('tight')

% topographic distribution
cfg                 = [];
cfg.layout          = 'neuromag306cmb.lay';
cfg.xlim            = timewin;
cfg.style           = 'straight';
cfg.comment         = 'no';
cfg.marker          = 'off';
cfg.colorbar        = 'southoutside';

cfg.figure = subplot(2,3,1);
ft_topoplotER(cfg, avg_loose_cmb{1});
title(['loose: ',conditions{1}]);

cfg.figure = subplot(2,3,2);
ft_topoplotER(cfg, avg_loose_cmb{2});
title(['loose: ',conditions{2}]);

 cfg.figure= subplot(2,3,3);
ft_topoplotER(cfg, avg_loose_cmb{3});
title(['loose: ',conditions{3}]);

cfg.figure = subplot(2,3,4);
ft_topoplotER(cfg, avg_tight_cmb{3});
title(['tight: ',conditions{1}]);

cfg.figure = subplot(2,3,5);
ft_topoplotER(cfg, avg_tight_cmb{2});
title(['tight: ',conditions{2}]);

 cfg.figure= subplot(2,3,6);
ft_topoplotER(cfg, avg_tight_cmb{3});
title(['tight: ',conditions{3}]);


%% Clean-up
rmpath('..')

