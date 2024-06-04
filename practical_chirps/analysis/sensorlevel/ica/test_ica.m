close all; clear all; clc;

% Settings
%--------------------------------------------------------------------------
% choose subjects
subjectlist = {'subject04'};

% choose files
% files2preproc = {'stories_maxfilter','olsa_maxfilter'};
files2preproc = 'all_maxfilter';
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);
% addpath for ica functions
addpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);

%% initialize
  
% subject selected out of codelist 
subject = subjectlist{1};
eval(subject)

filenames = get_filenames(subjectdata,files2preproc);
   
%% load and filter continuous data
cfg              = [];
cfg.dataset      = [subjectdata.rawdatadir filesep filenames{2} '.fif'];
% channels of interest and exclude bad channels
%cfg.channel        = [{'meg'},strcat('-',subjectdata.badchannels)]; 
cfg.channel      = 'meg';
cfg.demean       = 'yes';
cfg.detrend      = 'yes';
cfg.continuous   = 'yes';
cfg.coilaccuracy = 0;
cfg.bpfilter     = 'yes';
cfg.bpfreq       = [0.1,150];
cfg.dftfilter    = 'yes';        % enable notch filtering to eliminate power line noise
cfg.dftfreq      = [50 100 150]; % set up the frequencies for notch filtering
data_filtered    = ft_preprocessing(cfg);      
  
%% downsample data to reduce memory load
cfg              = [];
cfg.resamplefs   = 400;
cfg.detrend      = 'no';
data_filtered    = ft_resampledata(cfg,data_filtered);
    
%% ica on different sensortypes

% perform ica
cfg            = [];
% perform pca for maxfiltered meg data
%cfg.channel    = [{'meg'},strcat('-',subjectdata.badchannels)];
cfg.channel     = 'meg';
cfg.runica.pca = 64; % see natmeg tutorial for combined meg/eeg
cfg.method     = 'runica';
cfg.channel    = 'meg';
components     = ft_componentanalysis(cfg,data_filtered);  

%components.grad.chanunitold = data_filtered.grad.chanunit;
%components.grad.chantypeold = data_filtered.grad.chantype;
%components.grad.chanunit = data_filtered.grad.chanunit;
%components.grad.chantype = data_filtered.grad.chantype;

cfg           = [];
cfg.component = [1 2 3];
data_clean    = ft_rejectcomponent(cfg,components,data_filtered);

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])
