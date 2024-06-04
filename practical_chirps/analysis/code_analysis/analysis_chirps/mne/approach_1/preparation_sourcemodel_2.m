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

% resolution - specify only one hemisphere
resolution = '.L.midthickness.4k_fs_LR.surf.gii';
%resolution = '.L.midthickness.8k_fs_LR.surf.gii';

%--------------------------------------------------------------------------

% addpath for subject_files information
old = 'Z:'; new = '/media/till/Samsung_T5';

addpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new));
eval(subject)

%% Source model: Co-registration of the source space to the sensor-based head coordinate system

filename = [subjectdata.subjectname,resolution];
pathi    = replace([subjectdata.chirps_mne filesep 'sourcemodel' filesep],old,new);
path     = [pathi subjectdata.subjectname '_workbench' filesep filename];

% add right hemisphere
sourcemodel = ft_read_headshape({path, strrep(path,'.L.','.R.')});
sourcemodel = ft_determine_units(sourcemodel);

% load transformation matrices
transform_vox2acpc     = importdata([pathi,subjectdata.subjectname '_transform_vox2acpc.mat']);
transform_vox2neuromag = importdata([pathi,subjectdata.subjectname '_transform_vox2neuromag.mat']);

% concatenate transformation matrices (acpc -> neuromag)
transform_acpc2neuromag = transform_vox2neuromag/transform_vox2acpc;

% apply transformation
sourcemodel = ft_transform_geometry(transform_acpc2neuromag, sourcemodel);
%sourcemodel.inside = sourcemodel.atlasroi>0;

figure
ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', ...
             'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);
         
save([pathi,subjectdata.subjectname,'_sourcemodel.mat'],'sourcemodel');

%% Clean up
rmpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new))

