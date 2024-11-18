%--------------------------------------------------------------------------
% Plot sensorlevel results
%
% Check out: https://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging/
%--------------------------------------------------------------------------

close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

%% Script settings
%--------------------------------------------------------------------------
subject    = 'sub-01'; % 'sub-02'
conditions = settings.conditions;

%% Load data
%--------------------------------------------------------------------------
data     = importdata(fullfile(settings.path2project,'derivatives',subject,'sensorlevel',[subject,'_erf-N19mP30m.mat']));    
avg      = data.avg;
N_trials = data.N_trials;
clear data

%% Plot of event related fields for all sensors arranged topographically 
%--------------------------------------------------------------------------
% plot all 3 conditions together

% choose layout to look at gradiometers or magnetometers
% layout = 'neuromag306mag.lay';
layout = 'neuromag306planar.lay';

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbk';
ft_multiplotER(cfg, avg{:});
legend({conditions{1};conditions{2};conditions{3}});
title(subject)

%% Plot one specific channel
%--------------------------------------------------------------------------
% plot all 3 conditions together
chan2plot = 'MEG1221';

cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg{:});
legend({conditions{1};conditions{2};conditions{3}});
title(subject)

%% Plot the topographic distribution of the data
%--------------------------------------------------------------------------
timewin = [0.02 0.05]; % time window to look at in more detail

layout = 'neuromag306mag.lay';
% layout = 'neuromag306planar.lay';

cfg          = [];
cfg.layout   = layout;
cfg.xlim     = timewin;
cfg.style    = 'straight';
cfg.comment  = 'no';
cfg.marker   = 'off';
cfg.colorbar = 'southoutside';

cfg.figure = subplot(1,3,1);
ft_topoplotER(cfg, avg{1});
title(conditions{1});

cfg.figure = subplot(1,3,2);
ft_topoplotER(cfg, avg{2});
title(conditions{2});

cfg.figure= subplot(1,3,3);
ft_topoplotER(cfg, avg{3});
title(conditions{3});

%% Combine planar gradients of averaged data
%--------------------------------------------------------------------------
avg_cmb = cell(1,3);

for cidx=1:3
    cfg           = [];
    avg_cmb{cidx} = ft_combineplanar(cfg, avg{cidx});
end

cfg          = [];
cfg.layout   = 'neuromag306cmb.lay';
cfg.xlim     = timewin;
cfg.style    = 'straight';
cfg.comment  = 'no';
cfg.marker   = 'off';
cfg.colorbar = 'southoutside';

cfg.figure = subplot(1,3,1);
ft_topoplotER(cfg, avg_cmb{1});
title(conditions{1});

cfg.figure = subplot(1,3,2);
ft_topoplotER(cfg, avg_cmb{2});
title(conditions{2});

cfg.figure= subplot(1,3,3);
ft_topoplotER(cfg, avg_cmb{3});
title(conditions{3});

%% Clean-up
rmpath('..')