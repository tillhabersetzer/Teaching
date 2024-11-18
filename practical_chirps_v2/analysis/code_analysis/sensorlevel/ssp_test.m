%--------------------------------------------------------------------------
% Computation of event related fields over subjects (1,2,3,4) and 
% conditions {'clicks','upchirps','downchirps'}
%
% Check out: https://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging/
%
% Up to now:
% - ignore artifact channels (EOG, ECG), no ICA, no SSP applied
% - use simple z-value threshold for rejection of bad epochs, might be too 
%   easy and naiv
% - no headposition transformation and movement correction applied yet,
%   only simple tSSS method. The first point in particular may lead to 
%   questionable results when comparing different conditions because the
%   conditions were recorded in separate runs and hence have different head
%   positions.
%--------------------------------------------------------------------------

close all
clear 
clc 

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','subjectdata'))
eval('main_settings')

addpath(fullfile('..','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
subidx       = 3;
stimtype     = 'clicks';
TrigID.click = 1;
eventvalue   = TrigID.click;

subject  = ['sub-',num2str(subidx,'%02d')];
suffix   = '-raw_tsss';
% maxfilter discards ssp vectors
% datapath = fullfile(settings.path2bids,'derivatives',subject,'tsss',[subject,'_task-',stimtype,'_meg',suffix,'.fif']);
datapath = fullfile(settings.path2bids,subject,'meg',[subject,'_task-',stimtype,'_meg','.fif']);

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

% Demean trials
%--------------
cfg                = [];
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.05 0];
data               = ft_preprocessing(cfg,data);  

% Timelockanalysis
%-----------------
cfg = [];
avg = ft_timelockanalysis(cfg, data);

 
%% SSP

cfg          = [];
cfg.ssp      = 'all';
denoised_avg = ft_denoise_ssp(cfg,avg);

%% Comparison

cfg        = [];
cfg.layout = 'neuromag306mag_helmet.mat';
ft_multiplotER(cfg,avg, denoised_avg); 
       
