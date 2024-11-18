close all
clear 
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','..','subjectdata'))
eval('main_settings')

%% Script settings

% choose subject number
%--------------------------------------------------------------------------
subidx  = 4;
runidx  = 2;
subject = ['sub-',num2str(subidx,'%02d')];

% choose dataset for headshape extraction
megfile  = fullfile(settings.path2bids,subject,'meg',[subject,'_task-clicks_meg.fif']); 
mrifile  = fullfile(settings.path2bids,subject,'anat',[subject,'_run-0',num2str(runidx),'_T1w.nii']);
dir2save = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel');
%--------------------------------------------------------------------------

% make folder for data
%---------------------
if ~exist(dir2save, 'dir')
   mkdir(dir2save)
end

%% Preparation for Sourcemodel based on Surface Description
%--------------------------------------------------------------------------

% Extract Headshape
%------------------
headshape = ft_read_headshape(megfile,'unit','mm');

% Read in Mri
%------------
mri_orig = ft_read_mri(mrifile);

% 1. Coregistration to Neuromag Coordsystem based on Anatomical Landmarks
%------------------------------------------------------------------------
cfg            = [];
cfg.method     = 'interactive';
cfg.coordsys   = 'neuromag';
mri_realigned1 = ft_volumerealign(cfg, mri_orig);

% Incorporate Headshape Information (optional but recommended)
%-------------------------------------------------------------
cfg                       = [];
cfg.method                = 'headshape';
cfg.headshape.headshape   = headshape;
cfg.coordsys              = 'neuromag';
cfg.headshape.interactive = 'no';  % if manually: yes
cfg.headshape.icp         = 'yes'; % if manually: no / without automatic icp alignment - use only headshape
mri_realigned2            = ft_volumerealign(cfg, mri_realigned1);

% Check Results interactively in case of manual correction
%---------------------------------------------------------
cfg                       = [];
cfg.method                = 'headshape';
cfg.headshape.headshape   = headshape;
cfg.coordsys              = 'neuromag';
cfg.headshape.interactive = 'yes';
cfg.headshape.icp         = 'no'; 
mri_realigned3            = ft_volumerealign(cfg, mri_realigned2);

save(fullfile(dir2save,[subject,'_run-0',num2str(runidx),'_T1w-coregistered.mat']),'mri_realigned3'); % in mm


