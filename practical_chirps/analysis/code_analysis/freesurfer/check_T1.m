close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','subjectdata'))
eval('main_settings')

%% Define subjects
%--------------------------------------------------------------------------
subjects = 4;

% Load freesurfer-processed mri's
%--------------------------------

for subidx=subjects
     subject = ['sub-',num2str(subidx,'%02d')];
     mri_processed.(matlab.lang.makeValidName(subject)) = ft_read_mri(fullfile(settings.path2bids,'derivatives',subject,'freesurfer',[subject,'.nii']));
end

% Visualize data
%---------------
cfg = [];
for subidx=subjects
    subject = ['sub_',num2str(subidx,'%02d')];
    ft_sourceplot(cfg,mri_processed.(subject))
end

% Load original mri's
%--------------------

% only  first run
runidx = 1;

for subidx=subjects
     subject = ['sub-',num2str(subidx,'%02d')];
     filename = strcat(subject,'_run-0',num2str(runidx),'_T1w.nii');  
     mri_orig.(matlab.lang.makeValidName(subject)) = ft_read_mri(fullfile(settings.path2bids,subject,'anat',filename));
end

% Visualize data
%---------------
cfg = [];
for subidx=subjects
    subject = ['sub_',num2str(subidx,'%02d')];
    ft_sourceplot(cfg,mri_orig.(subject))
end