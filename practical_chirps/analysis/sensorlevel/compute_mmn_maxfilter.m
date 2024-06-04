close all
clear 
clc 

% checkout: https://www.fieldtriptoolbox.org/workshop/natmeg/preprocessing/

% Results
% Application of maxfilter results in significant differences only for
% magnetometers. In that case maxfilter helps to reduce artifacts.
% The gradiometers look almost the same.

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

%% Script settings
%--------------------------------------------------------------------------
subidx = [2];
subject  = ['sub-',num2str(subidx,'%02d')];

%--------------------------------------------------------------------------

% Path for data
%--------------------------------------------------------------------------
datapath{1} = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-mmn.fif']);
datapath{2} = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-mmn-raw_tsss.fif']);

% data overview
%--------------------------------------------------------------------------
cfg            = [];
cfg.dataset    = datapath{2};
cfg.continuous = 'yes';
cfg.channel    = 'megmag';
cfg.viewmode   = 'vertical';
cfg.blocksize  = 1;                             
ft_databrowser(cfg);

cfg            = [];
cfg.dataset    = datapath{2};
cfg.continuous = 'yes';
cfg.channel    = 'megplanar';
cfg.viewmode   = 'vertical';
cfg.blocksize  = 1;                             
ft_databrowser(cfg);

cfg                = [];
cfg.dataset        = datapath{2};
cfg.continuous     = 'yes';
cfg.channel        = 'EEG';
cfg.viewmode       = 'vertical';
cfg.blocksize      = 1;   
cfg.preproc.demean = 'yes';    
ft_databrowser(cfg);

%% Computation of evoked fields
%--------------------------------------------------------------------------

avg       = cell(2,3);
trialinfo = zeros(2,2);

for cidx = 1:2 % without/with maxfilter

    cfg              = [];
    cfg.dataset      = datapath{cidx};
    cfg.channel      = {'meg'}; 
    cfg.continuous   = 'yes';
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = [1,30];
    cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
    data             = ft_preprocessing(cfg); 

    % Define trials
    %--------------
    cfg                     = [];
    cfg.dataset             = datapath{cidx};
    cfg.trialfun            = 'ft_trialfun_general'; 
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = [1,2]; % standard / deviant
    cfg.trialdef.prestim    = 0.2;                  
    cfg.trialdef.poststim   = 0.8;                 
    cfg                     = ft_definetrial(cfg);
    trl                     = cfg.trl;
 
    % Epoch data
    %-----------
    cfg     = [];
    cfg.trl = trl;  
    epoched = ft_redefinetrial(cfg,data); 

    % Demean trials
    %--------------
    cfg                = [];
    cfg.demean         = 'yes';
    cfg.baselinewindow = [-0.2 0];
    epoched            = ft_preprocessing(cfg,epoched);  

    % Semi automatic artifact rejection
    %----------------------------------
    % separately for magnetometers
    cfg             = [];
    cfg.metric      = 'zvalue';
    cfg.layout      = 'neuromag306all.lay';
    cfg.channel     = 'megmag';
    cfg.keepchannel = 'yes';  % This keeps those channels that are not displayed in the data
    epoched         = ft_rejectvisual(cfg,epoched);
    
    % separately for gradiometers
    cfg.channel = 'megplanar';
    epoched     = ft_rejectvisual(cfg,epoched);

    % Downsample data
    %----------------
    cfg            = [];
    cfg.resamplefs = 100;
    epoched        = ft_resampledata(cfg,epoched); 

    % Timelockanalysis
    %-----------------
    % standard tones
    cfg         = [];
    cfg.trials  = find(epoched.trialinfo(:,1) == 1);
    avg{cidx,1} = ft_timelockanalysis(cfg, epoched);

    % deviant tones
    cfg         = [];
    cfg.trials  = find(epoched.trialinfo(:,1) == 2);
    avg{cidx,2} = ft_timelockanalysis(cfg, epoched);

    clear epoched
    trialinfo(cidx,1) = length(avg{cidx,1}.cfg.trials);
    trialinfo(cidx,2) = length(avg{cidx,2}.cfg.trials);

    % differences
    cfg           = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    cfg.operation = 'x1 - x2'; % see https://www.fieldtriptoolbox.org/workshop/natmeg/dipolefitting/
    avg{cidx,3}   = ft_math(cfg, avg{cidx,2}, avg{cidx,1});
  
end

%% Visualize data - timeseries

 % gradiometers look way superior
 
layout = 'neuromag306mag.lay';
% layout = 'neuromag306planar.lay';

% without maxfilter
%------------------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{1,1}, avg{1,2}, avg{1,3});
legend({'standard','deviant','difference'});
subtitle('without maxfilter')

% with maxfilter
%---------------
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{2,1}, avg{2,2}, avg{2,3});
legend({'standard','deviant','difference'});
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
subtitle('standard');

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{1,2}, avg{2,2});
legend({'without maxfilter','with maxfilter'})
subtitle('deviant');

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'br';
ft_multiplotER(cfg, avg{1,3}, avg{2,3});
legend({'without maxfilter','with maxfilter'})
subtitle('difference');

%% Visualize data - topography

layout = 'neuromag306mag.lay';
% layout = 'neuromag306planar.lay';

cfg                 = [];
cfg.layout          = layout;
cfg.xlim            = [0.08 0.15];
cfg.style           = 'straight';
cfg.comment         = 'no';
cfg.marker          = 'off';
cfg.colorbar        = 'southoutside';

cfg.figure = subplot(2,3,1);
ft_topoplotER(cfg, avg{1,1});
title('without mf: standard');

cfg.figure = subplot(2,3,2);
ft_topoplotER(cfg, avg{1,2});
title('without mf: deviant');

 cfg.figure= subplot(2,3,3);
ft_topoplotER(cfg, avg{1,3});
title('without mf: difference');

cfg.figure = subplot(2,3,4);
ft_topoplotER(cfg, avg{2,1});
title('with mf: standard');

cfg.figure = subplot(2,3,5);
ft_topoplotER(cfg, avg{2,2});
title('with mf: deviant');

 cfg.figure= subplot(2,3,6);
ft_topoplotER(cfg, avg{2,3});
title('with mf: difference');

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
legend({'standard','deviant','difference'});
subtitle('without maxfilter')

cfg          = [];
cfg.layout   = 'neuromag306cmb.lay';
ft_multiplotER(cfg, avg_cmb{2,1}, avg_cmb{2,2}, avg_cmb{2,3});
legend({'standard','deviant','difference'});
subtitle('with maxfilter')

% topographic distribution
cfg                 = [];
cfg.layout          = 'neuromag306cmb.lay';
cfg.xlim            = [0.08 0.15];
cfg.style           = 'straight';
cfg.comment         = 'no';
cfg.marker          = 'off';
cfg.colorbar        = 'southoutside';

cfg.figure = subplot(2,3,1);
ft_topoplotER(cfg, avg_cmb{1,1});
title('without mf: standard');

cfg.figure = subplot(2,3,2);
ft_topoplotER(cfg, avg_cmb{1,2});
title('without mf: deviant');

 cfg.figure= subplot(2,3,3);
ft_topoplotER(cfg, avg_cmb{1,3});
title('without mf: difference');

cfg.figure = subplot(2,3,4);
ft_topoplotER(cfg, avg_cmb{2,1});
title('with mf: standard');

cfg.figure = subplot(2,3,5);
ft_topoplotER(cfg, avg_cmb{2,2});
title('with mf: deviant');

 cfg.figure= subplot(2,3,6);
ft_topoplotER(cfg, avg_cmb{2,3});
title('with mf: difference');

%%
rmpath('..')







