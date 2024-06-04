close all
clear 
clc 

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

%% Script settings
%--------------------------------------------------------------------------
subidx = [2];

conditions = settings.conditions;
%--------------------------------------------------------------------------

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

%% Preprocess data
%--------------------------------------------------------------------------

avg      = cell(2,3);
N_trials = zeros(2,3);

for cidx = 1:3 % loop over conditions

    % Path for MEG data
    %------------------
    subject   = ['sub-',num2str(subidx,'%02d')];
    condition = conditions{cidx};

    datapath{1} = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);
    datapath{2} = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-',condition,'-raw_tsss.fif']);

    switch condition
        case 'clicks'
            eventvalue = TrigID.click;
        case 'upchirps'
            eventvalue = TrigID.upchirp;
        case 'downchirps'
            eventvalue = TrigID.downchirp;
    end

    for midx = 1:2 % loop over maxfilter

        % Define trials
        %--------------
        cfg                     = [];
        cfg.dataset             = datapath{midx};
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
        cfg.dataset      = datapath{midx};
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

        % Demean trials
        %--------------
        cfg                = [];
        cfg.demean         = 'yes';
        cfg.baselinewindow = [-0.05 0];
        epoched            = ft_preprocessing(cfg,epoched);  
    
        % Semi automatic artifact rejection
        %----------------------------------
        % separately for magnetometers
%             cfg             = [];
%             cfg.metric      = 'zvalue';
%             cfg.layout      = 'neuromag306all.lay';
%             cfg.channel     = 'megmag';
%             cfg.keepchannel = 'yes';  % This keeps those channels that are not displayed in the data
%             epoched         = ft_rejectvisual(cfg,epoched);
        
        % separately for gradiometers
%             cfg.channel = 'megplanar';
%             epoched     = ft_rejectvisual(cfg,epoched);

        % Timelockanalysis
        %-----------------
        cfg            = [];
        avg{midx,cidx} = ft_timelockanalysis(cfg, epoched);

        % Note down amount of preserved epochs
        %-------------------------------------
        N_trials(midx,cidx) = length(avg{midx,cidx}.cfg.previous.trials);

        clear data epoched

    end

end

%% Plot of event related fields for all sensors arranged topographically 
%--------------------------------------------------------------------------
% layout = 'neuromag306mag.lay';
layout = 'neuromag306planar.lay';

% without maxfilter
%------------------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{1,1}, avg{1,2}, avg{1,3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('without maxfilter')

% with maxfilter
%---------------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{2,1}, avg{2,2}, avg{2,3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('with maxfilter')

% comparison
%-----------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{1,1}, avg{2,1});
legend({'without maxfilter','with maxfilter'})
subtitle(conditions{1});

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{1,2}, avg{2,2});
legend({'without maxfilter','with maxfilter'})
subtitle(conditions{2});

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{1,3}, avg{2,3});
legend({'without maxfilter','with maxfilter'})


%% Plot one specific channel
%--------------------------------------------------------------------------
% plot all 3 conditions together
chan2plot = 'MEG1411';

cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg{1,:});
legend({conditions{1};conditions{2};conditions{3}});
title(['without maxfilter: ',subject])

cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg{2,:});
legend({conditions{1};conditions{2};conditions{3}});
title(['with maxfilter: ',subject])

%% Visualize data - topography

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
ft_topoplotER(cfg, avg{1,1});
title(['without mf: ',conditions{1}]);

cfg.figure = subplot(2,3,2);
ft_topoplotER(cfg, avg{1,2});
title(['without mf: ',conditions{2}]);

 cfg.figure= subplot(2,3,3);
ft_topoplotER(cfg, avg{1,3});
title(['without mf: ',conditions{3}]);

cfg.figure = subplot(2,3,4);
ft_topoplotER(cfg, avg{2,1});
title(['with mf: ',conditions{1}]);

cfg.figure = subplot(2,3,5);
ft_topoplotER(cfg, avg{2,2});
title(['with mf: ',conditions{2}]);

 cfg.figure= subplot(2,3,6);
ft_topoplotER(cfg, avg{2,3});
title(['with mf: ',conditions{3}]);

%% Combining planar gradiometers
% To help with identifying underlying sources we should make use of the 
% other channels in the data. The planar gradiometers are often more easily 
% interpreted, because they are most sensitive right above a source. 

% Combine planar
avg_cmb = cell(2,3);
cfg     = [];
for n = 1:6
    avg_cmb{n} = ft_combineplanar(cfg, avg{n});
end

% Plot results
%-------------
cfg        = [];
cfg.layout = 'neuromag306cmb.lay';
ft_multiplotER(cfg, avg_cmb{1,1}, avg_cmb{1,2}, avg_cmb{1,3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('without maxfilter')

cfg          = [];
cfg.layout   = 'neuromag306cmb.lay';
ft_multiplotER(cfg, avg_cmb{2,1}, avg_cmb{2,2}, avg_cmb{2,3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('with maxfilter')

% topographic distribution
cfg                 = [];
cfg.layout          = 'neuromag306cmb.lay';
cfg.xlim            = timewin;
cfg.style           = 'straight';
cfg.comment         = 'no';
cfg.marker          = 'off';
cfg.colorbar        = 'southoutside';

cfg.figure = subplot(2,3,1);
ft_topoplotER(cfg, avg_cmb{1,1});
title(['without mf: ',conditions{1}]);

cfg.figure = subplot(2,3,2);
ft_topoplotER(cfg, avg_cmb{1,2});
title(['without mf: ',conditions{2}]);

cfg.figure= subplot(2,3,3);
ft_topoplotER(cfg, avg_cmb{1,3});
title(['with mf: ',conditions{3}]);

cfg.figure = subplot(2,3,4);
ft_topoplotER(cfg, avg_cmb{2,1});
title(['with mf: ',conditions{1}]);

cfg.figure = subplot(2,3,5);
ft_topoplotER(cfg, avg_cmb{2,2});
title(['with mf: ',conditions{2}]);

 cfg.figure= subplot(2,3,6);
ft_topoplotER(cfg, avg_cmb{2,3});
title(['with mf: ',conditions{3}]);


%% Cleaning
rmpath('..')
