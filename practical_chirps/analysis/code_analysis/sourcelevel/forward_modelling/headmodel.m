%--------------------------------------------------------------------------
% Within this script a singleshell headmodel for MEG is calculated. In
% addition the mri is preprocessed for the freesurfer pipeline that creates
% surface based sourcemodels.
%--------------------------------------------------------------------------

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
subject = ['sub-',num2str(subidx,'%02d')];

% choose dataset for headshape extraction
megfile  = fullfile(settings.path2bids,subject,'meg',[subject,'_task-clicks_meg.fif']); 
% choose mri file
% mrifile  = fullfile(settings.path2bids,subject,'anat',[subject,'_T1w.nii']);
mrifile  = fullfile(settings.path2bids,'derivatives',subject,'freesurfer',[subject,'.nii']);

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
mri_coreg                 = ft_volumerealign(cfg, mri_realigned2);

% -> freesurfer mri's have already been resliced

% 2. Reslice
%-----------
% Check if full head is still in mri after reslicing, might fail for big
% heads and is important in case you may need a bem headmodel which
% includes the scalp
% cfg            = [];
% cfg.resolution = 1;
% cfg.dim        = [256 256 256];
% cfg.yrange     = [-130,125]; % could be adjusted, see above remark, sum up to 255 -> 256 voxel
% cfg.coordsys   = 'neuromag';
% mri_resliced   = ft_volumereslice(cfg,mri_realigned3); 

% Check 
%------
cfg = [];
ft_sourceplot(cfg, mri_coreg);
% cfg = [];
% ft_sourceplot(cfg, mri_resliced);

% Save the Resliced Anatomy in a FreeSurfer Compatible Format
%------------------------------------------------------------
% Volumewrite saves no Information about the Coordinate System!
% cfg           = [];
% cfg.filename  = fullfile(dir2save,[subject,'_T1w-coregistered.mgz']);
% cfg.filetype  = 'mgz';
% cfg.parameter = 'anatomy';
% ft_volumewrite(cfg, mri_resliced);
% --> Needed for Headmodel and FreeSurfer

% 3. Segmentation (most time consuming)
%--------------------------------------
cfg           = [];
cfg.output    = {'brain','skull','scalp'}; % you only need brain for meg
% mri_segmented = ft_volumesegment(cfg, mri_resliced);
mri_segmented = ft_volumesegment(cfg, mri_coreg);

% Copy the Anatomy into the segmented Mri
% mri_segmented.anatomy = mri_resliced.anatomy;
mri_segmented.anatomy = mri_coreg.anatomy;

% Check
%------
cfg              = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, mri_segmented);

% 4. Headmodel
%-------------
cfg        = [];
cfg.method = 'singleshell';
vol        = ft_prepare_headmodel(cfg, mri_segmented);

% Save data
%----------
save(fullfile(dir2save,[subject,'_T1w-coregistered.mat']),'mri_coreg'); % in mm
% save(fullfile(dir2save,[subject,'_T1w-resliced.mat']),'mri_resliced'); % in mm
save(fullfile(dir2save,[subject,'_T1w-segmented.mat']),'mri_segmented'); % in mm
save(fullfile(dir2save,[subject,'_headmodel-singleshell.mat']),'vol'); % in mm

% Save coregistration
%--------------------
addpath(fullfile('..','..','helper_functions'))

cfg           = [];
cfg.shape     = shape;
cfg.transform = mri_coreg.transform;
[pos, label]  = landmarks_hcs2neuromag(cfg);

landmarks.lpa = pos(strcmp(label,'LPA'),:);
landmarks.rpa = pos(strcmp(label,'RPA'),:);
landmarks.nas = pos(strcmp(label,'Nasion'),:);
% landmarks.zpoint = [120,120,250];

save(fullfile(dir2save,[subject,'_landmarks-coregistered.mat']),'landmarks'); % in mm

% Check
%------
% Same Units for Plot
grad   = ft_read_sens(megfile,'senstype','meg'); % cm
shape  = ft_read_headshape(megfile,'unit','cm');
vol_cm = ft_convert_units(vol,'cm');

figure
ft_plot_headmodel(vol_cm);
hold on
ft_plot_sens(grad, 'style', '*b');
ft_plot_headshape(shape);

