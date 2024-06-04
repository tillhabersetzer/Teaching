close all; clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!SCRIPT IS DESIGNED FOR LINUX!!!!
% to work on windows replace new = '/media/till/Samsung_T5' with new = old
% then everthing with the pays stays the same
%
% and you have to change mgz to nifti format on windows!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subject  = 'subject04';

% choose dataset for headshape extraction
filename = 'fixer_olsa_tsss_mc'; 
%--------------------------------------------------------------------------

% addpath for subject_files information
old = 'Z:'; new = '/media/till/Samsung_T5';

addpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new));
eval(subject)

% make folder for data
new_dir = replace([subjectdata.chirps_mne filesep 'sourcemodel'],old,new);
if ~exist(new_dir, 'dir')
   mkdir(new_dir)
end

%% preparation for sorucemodel based on surface description
%----------------------------------------------------------

% extract headshape
%------------------
dataset = replace([subjectdata.rawdatadir filesep filename '.fif'],old,new);
shape   = ft_read_headshape(dataset,'unit','m');

% read in mri
%------------
mri_orig = ft_read_mri(replace(subjectdata.mri,old,new));

% coregistration to neuromag coordsys
%------------------------------------
cfg            = [];
cfg.method     = 'interactive';
cfg.coordsys   = 'neuromag';
mri_realigned1 = ft_volumerealign(cfg, mri_orig);

cfg            = [];
cfg.method     = 'headshape';
cfg.headshape  = shape;
cfg.coordsys   = 'neuromag';
mri_realigned2 = ft_volumerealign(cfg, mri_realigned1);

% reslice
%-------- 
cfg            = [];
cfg.resolution = 1;
cfg.dim        = [256 256 256];
cfg.coordsys   = 'neuromag';
mri_resliced   = ft_volumereslice(cfg,mri_realigned2); 

% check 
%------
cfg = [];
ft_sourceplot(cfg, mri_realigned2);
cfg = [];
ft_sourceplot(cfg, mri_resliced);

transform_vox2neuromag = mri_resliced.transform;
save(replace(fullfile(new_dir,sprintf('%s_transform_vox2neuromag',subjectdata.subjectname)),old,new),'transform_vox2neuromag')

% save the resliced anatomy in a FreeSurfer compatible format
%------------------------------------------------------------
% volumewrite saves no information about the coordinate system!
cfg           = [];
cfg.filename  = replace(fullfile(new_dir,sprintf('%s_neuromag.mgz',subjectdata.subjectname)),old,new);
cfg.filetype  = 'mgz';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_resliced);
% --> needed for headmodel

% coregistration to acpc coordsys
%--------------------------------
cfg          = [];
cfg.method   = 'interactive';
cfg.coordsys = 'acpc';
mri_acpc     = ft_volumerealign(cfg, mri_resliced);

transform_vox2acpc = mri_acpc.transform;
save(replace(fullfile(new_dir,sprintf('%s_transform_vox2acpc',subjectdata.subjectname)),old,new),'transform_vox2acpc')

% save the resliced anatomy in a FreeSurfer compatible format
%------------------------------------------------------------
cfg           = [];
cfg.filename  = replace(fullfile(new_dir,sprintf('%s.mgz',subjectdata.subjectname)),old,new);
cfg.filetype  = 'mgz';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_acpc);
% --> needed for freesurfer

%% Clean up
rmpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new))
