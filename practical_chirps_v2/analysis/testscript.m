close all
clear all
clc

% Filter continuous data to avoid edge artifacts
%-----------------------------------------------
cfg          = [];
cfg.dataset  = 'M:\Blockpraktikum2022\results\sub-01\meg\upchirps.fif';
cfg.channel  = 'meg';
cfg.bpfilter = 'yes';
cfg.bpfreq   = [16,120];
data         = ft_preprocessing(cfg); 

% figure
% plot(data.time{1},data.trial{1}(382,:))

% Define trials
%--------------
cfg                     = [];
cfg.dataset             = 'M:\Blockpraktikum2022\results\sub-01\meg\upchirps.fif';
cfg.trialfun            = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype  = 'STI101';
cfg.trialdef.eventvalue = 2; 
cfg.trialdef.prestim    = 0.05;                  % in seconds
cfg.trialdef.poststim   = 0.35;                  % in seconds
cfg                     = ft_definetrial(cfg);
trl                     = cfg.trl;

% Epoch data
%-----------
cfg          = [];
cfg.trl      = trl;  
data_epoched = ft_redefinetrial(cfg,data); 

% % Demean trials
% %--------------
% cfg                = [];
% cfg.demean         = 'yes';
% cfg.baselinewindow = [-0.05 0];
% data               = ft_preprocessing(cfg,data);  

cfg = [];
avg = ft_timelockanalysis(cfg, data_epoched);

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
ft_multiplotER(cfg, avg);


shape = ft_read_headshape('M:\Blockpraktikum2022\results\sub-01\meg\upchirps.fif');
ft_plot_headshape(shape)







