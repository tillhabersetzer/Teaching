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
%--------------------------------------------------------------------------

% addpath for subject_files information
old = 'Z:'; new = '/media/till/Samsung_T5';

addpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new));
eval(subject)

%% generate headmodel

pathi = replace([subjectdata.chirps_mne filesep 'sourcemodel' filesep],old,new);
path  = [pathi,sprintf('%s_neuromag.mgz',subjectdata.subjectname)];

% load already aligned and resliced mri from preparation_sourcemodel_1
mri          = ft_read_mri(path);
mri.coordsys = 'neuromag';

% segmentation (most time consuming)
%-----------------------------------
cfg           = [];
cfg.output    = {'brain'}; % you only need brain for meg
mri_segmented = ft_volumesegment(cfg, mri);

% copy the anatomy into the segmented mri
mri_segmented.anatomy = mri.anatomy;

% check
%------
cfg              = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, mri_segmented);

% notice
% you can put the mri_segmented directly into the ft_prepare_headmodel
% function, but it is also possible to put a mesh grid into there

% triangulated meshes
%--------------------
cfg             = [];
cfg.tissue      = 'brain';
cfg.method      = 'projectmesh';
cfg.numvertices = 3000;
% cfg.method      = 'isosurface'; 
% cfg.numvertices = 16000;
mesh_brain      = ft_prepare_mesh(cfg, mri_segmented);

% check
%------
figure
ft_plot_mesh(mesh_brain, 'edgecolor', 'none', 'facecolor', 'skin')
material dull
camlight
lighting phong

% headmodel
%----------
cfg        = [];
cfg.method = 'singleshell';
vol        = ft_prepare_headmodel(cfg, mesh_brain);
% vol        = ft_prepare_headmodel(cfg, mri_segmented);

save([pathi,subjectdata.subjectname,'_headmodel_meg.mat'],'vol'); % in mm

% check
%------
dataset = replace([subjectdata.rawdatadir filesep 'fixer_olsa_tsss_mc.fif'],old,new);
% same units for plotting
vol_cm  = ft_convert_units(vol,'cm');
grad    = ft_read_sens(dataset,'senstype','meg'); % cm
shape   = ft_read_headshape(dataset,'unit','cm');

figure
ft_plot_headmodel(vol_cm);
hold on
ft_plot_sens(grad, 'style', '*b');
ft_plot_headshape(shape);

%% Clean up
rmpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new))
