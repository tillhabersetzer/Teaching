close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the avarages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose only one file
vp        = 'ck06an19_vp';
%file      = 'fixer_olsa_tsss_mc';
file      = 'story1part2_tsss_mc';
path_root = ['Z:' filesep 'data_megproject' filesep vp filesep];
%--------------------------------------------------------------------------

path_dataset      = [path_root filesep 'recordings' filesep file '.fif'];
addpath(['Z:' filesep 'validation_pilot_update' filesep 'helper_functions']);

%event = ft_read_event(path_dataset{n}, 'checkmaxfilter', 'no');


%% define trials
cfg                     = [];
cfg.dataset             = path_dataset;
cfg.trialfun            = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype  = 'STI101';
cfg.trialdef.eventvalue = 2; % 
cfg.trialdef.prestim    = 0.05; % in seconds
cfg.trialdef.poststim   = 0.35; % in seconds
cfg                     = ft_definetrial(cfg);
trl                     = cfg.trl;

%% detect artifacts
    cfg                       = [];
    cfg.trl                   = trl; 
    cfg.artfctdef.eog.channel = {'eog'};
    cfg.dataset               = path_dataset;
    [~, artifact_eog]         = ft_artifact_eog(cfg);
    
    cfg                = [];
    cfg.trl            = trl;
    cfg.dataset        = path_dataset;
    [~, artifact_jump] = ft_artifact_jump(cfg);
    
%     cfg                  = [];
%     cfg.trl              = trl;
%     cfg.dataset          = path_dataset;
%     [~, artifact_muscle] = ft_artifact_muscle(cfg); 
    
    cfg                              = [];
    cfg.trl                          = trl;
    cfg.dataset                      = path_dataset;
    cfg.artfctdef.eog.artifact       = artifact_eog;
    cfg.artfctdef.jump.artifact      = artifact_jump;
    % cfg.artfctdef.muscle.artifact    = artifact_muscle;
    cfg = ft_rejectartifact(cfg);
    trl_new = cfg.trl;

%     
%     cfg                = [];
%     cfg.trl            = trl;
%     cfg.dataset        = path_dataset{n};
%     [~, artifact_jump] = ft_artifact_jump(cfg);
%     
%     cfg                  = [];
%     cfg.trl              = trl;
%     cfg.dataset          = path_dataset{n};
%     [~, artifact_muscle] = ft_artifact_muscle(cfg); 

%       cfg                              = [];
%       cfg.trl                          = trl;
%       cfg.dataset                      = path_dataset{n};
%       cfg.artfctdef.threshold.bpfilter = 'no'; 
%       cfg.artfctdef.threshold.channel  = {'megmag', '-MEG2143'};
%       cfg.artfctdef.threshold.min      = -1000*10^-15; % T
%       cfg.artfctdef.threshold.max      =  1000*10^15;
%       cfg.artfctdef.threshold.range    = 3000*10^-15;
%       [~, artifact_thres_mag]          = ft_artifact_threshold(cfg); 
      
%       cfg.artfctdef.threshold.channel  = {'megplanar', '-MEG2143'};
%      cfg.artfctdef.threshold.min      = -1000*10^-15*10^-2; % T
%      cfg.artfctdef.threshold.max      =  1000*10^15*10^-2;
%       cfg.artfctdef.threshold.range    = 3000*10^-15/0.04;
%       [~, artifact_thres_planar]       = ft_artifact_threshold(cfg); 
      
%       artifact_thres = [artifact_thres_mag;artifact_thres_planar];
       
      % check
      % A = max(abs(data.trial{1}),[],2);

% 
%      cfg                            = [];
%      cfg.trl                        = trl;
%      cfg.dataset                    = path_dataset{n};
%      cfg.artfctdef.eog.artifact    = artifact_eog;
%      cfg.artfctdef.jump.artifact   = artifact_jump;
%      cfg.artfctdef.muscle.artifact = artifact_muscle;
%      cfg.artfctdef.threshold.artifact = artifact_thres;
%      cfg = ft_rejectartifact(cfg);
%      trl = cfg.trl;

%end

%% preprocess the data
cfg                = [];
cfg.trl            = trl_new;
cfg.dataset        = path_dataset;
cfg.channel        = {'MEG', '-MEG2143'};    
cfg.demean         = 'yes';
cfg.detrend        = 'yes';
cfg.baselinewindow = [-0.05 0];
cfg.bpfilter       = 'yes';
cfg.bpfreq         = [1,120];
cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
cfg.dftfreq        = [60 120 180]; % set up the frequencies for notch filtering
data               = ft_preprocessing(cfg);

%% average
cfg = [];
avg = ft_timelockanalysis(cfg,data);

%% plot
figure('name','mag')
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
cfg.xlim       = [0, 0.2];
ft_multiplotER(cfg, avg);

figure('name','grad')
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
ft_multiplotER(cfg, avg);

figure
cfg.channel = 'MEG1612';
ft_singleplotER(cfg,avg);

figure
cfg          = [];
cfg.xlim     = [0.1 0.3];
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306mag.lay';
ft_topoplotER(cfg,avg);

figure
cfg        = [];
cfg.xlim   = [0 : 0.1 : 0.4];  % Define 12 time intervals
cfg.layout = 'neuromag306mag.lay';
ft_topoplotER(cfg,avg);

%% planar gradient
cfg        = [];
cfg.method = 'sum';
avg_comb   = ft_combineplanar(cfg,avg);

figure
cfg        = [];
cfg.xlim   = [0 : 0.1 : 0.4];  % Define 12 time intervals
cfg.layout = 'neuromag306cmb.lay';
ft_topoplotER(cfg,avg_comb);

figure('name','cmb')
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg, avg_comb);




%% have a look
label = importdata(['Z:' filesep 'validation_pilot_update' filesep 'helper_functions' filesep 'meg_label.mat']);

cfg                        = [];
cfg.dataset                = path_dataset;
cfg.trl                    = trl;
cfg.plotevents             = 'no';
cfg.artfctdef.eog.artifact = artifact_eog;
desired = {'MF-M1'}; % frontral magnetometer
%desired = {'MF-MG'}; % frontral gradiometer
%desired = {'M1'};    % only Magnetometer
my_selection               = channel_selection(desired,label);
cfg.channel                = ['eog';'ecg';my_selection];
cfg.viewmode               = 'vertical';
cfg.blocksize              = 50;
ft_databrowser(cfg);

cfg                = [];
cfg.dataset        = path_dataset;
% channels of interest and exclude bad channels
cfg.channel        = {'eog','ecg','MEG', '-MEG2143'}; 
cfg.demean         = 'yes';
cfg.detrend        = 'yes';
cfg.continuous     = 'yes';
cfg.bpfilter       = 'yes';
cfg.bpfreq         = [0.1,200];
cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
cfg.dftfreq        = [60 120 180]; % set up the frequencies for notch filtering
data_filtered      = ft_preprocessing(cfg);

% cfg                       = [];
% cfg.artfctdef.eog.channel = {'eog'};
% cfg.dataset               = path_dataset;
% [~, artifact_eog2]         = ft_artifact_eog(cfg,data_filtered);


% cfg                        = [];
% %cfg.trl                    = trl;
% cfg.plotevents             = 'no';
% cfg.artfctdef.eog.artifact = artifact_eog;
% cfg.artfctdef.jump.artifact = artifact_eog2;
% desired = {'MF-M1'}; % frontral magnetometer
% %desired = {'MF-MG'}; % frontral gradiometer
% %desired = {'M1'};    % only Magnetometer
% my_selection               = channel_selection(desired,label);
% cfg.channel                = ['eog';'ecg';my_selection];
% cfg.viewmode               = 'vertical';
% cfg.blocksize              = 100;
% ft_databrowser(cfg,data_filtered);


% Artefaktdetektion von Augenartefakten hängt von der Anzahl der Daten und
% der darin vorhandenen Artefakte ab. Dies liegt am Grenzwert, der durch
% die z-Transformation bestimmt wird. Gibt es mehr große Artefakte liegt
% der Schwellenwert höher.


