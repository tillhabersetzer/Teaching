%--------------------------------------------------------------------------
% https://www.fieldtriptoolbox.org/example/coregistration_quality_control/
%--------------------------------------------------------------------------

close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','..'))
eval('main_settings')

% check whether maxfiltered data should be analyzed
maxfilter = settings.maxfilter;

%% Script settings

% choose subject number
%--------------------------------------------------------------------------
subidx  = 2;
subject = ['sub-',num2str(subidx,'%02d')];

% choose dataset for headshape extraction
if maxfilter
    megfile  = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-clicks-raw_tsss.fif']);
    megfile2 = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-upchirps-raw_tsss.fif']);
    megfile3 = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-downchirps-raw_tsss.fif']);
else
    megfile  = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-clicks.fif']); 
    megfile2 = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-upchirps.fif']); 
    megfile3 = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-downchirps.fif']); 
end

mrifile  = fullfile(settings.path2project,'rawdata',subject,'anat',[subject,'_T1w.nii']);

%--------------------------------------------------------------------------

%% Load required data
%--------------------------------------------------------------------------
headmodel_meg   = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_headmodel-singleshell.mat']));
mri             = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_T1w-segmented.mat']));
sourcemodel_meg = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_sourcemodel-volumetric.mat']));
headshape       = ft_read_headshape(megfile,'unit','mm');

% Same Units for Plot (mm)
grad  = ft_read_sens(megfile,'senstype','meg');  % cm
grad2 = ft_read_sens(megfile2,'senstype','meg'); % cm
grad3 = ft_read_sens(megfile3,'senstype','meg'); % cm
grad  = ft_convert_units(grad,'mm');
grad2 = ft_convert_units(grad2,'mm');
grad3 = ft_convert_units(grad3,'mm');

%% figure 1, All together

figure
ft_plot_sens(grad)
ft_plot_headshape(headshape)
ft_plot_headmodel(headmodel_meg)
ft_plot_axes([], 'unit', 'mm');

%% figure 2, MRI anatomy and brain segmentation

cfg              = [];
cfg.anaparameter = 'anatomy';
cfg.funparameter = 'brain';
% cfg.funparameter = 'scalp';
cfg.location     = [0 0 60];
ft_sourceplot(cfg, mri)

%% figure 3 and 4, MRI anatomy and headmodel

location = [0 0 60];
figure
ft_plot_ortho(mri.anatomy, 'transform', mri.transform, 'location', location, 'intersectmesh', headmodel_meg.bnd,'unit','mm')

figure
ft_plot_montage(mri.anatomy, 'transform', mri.transform, 'intersectmesh', headmodel_meg.bnd,'unit','mm')

%% figure 5, MRI scalp surface and polhemus headshape
% problems with surface extraction due to "worse" mri segmentation likely
cfg             = [];
cfg.tissue      = {'scalp'};
cfg.method      = 'isosurface';
cfg.numvertices = 10000;
scalp           = ft_prepare_mesh(cfg, mri);

figure
ft_plot_mesh(scalp, 'facecolor', 'skin')
lighting phong
camlight left
camlight right
material dull
alpha 0.5
ft_plot_headshape(headshape, 'vertexcolor', 'k');

%% figure 6, MRI and anatomical landmarks
% seems not to work :/
figure
for i=1:3
  subplot(2,2,i)
  figure
  title(headshape.fid.label{i});
  location = headshape.fid.pos(i,:);
  ft_plot_ortho(mri.anatomy, 'transform', mri.transform, 'style', 'intersect', 'location', location, 'plotmarker', location, 'markersize', 5, 'markercolor', 'y')
end

%% figure 7, MRI scalp surface and anatomical landmarks

figure
ft_plot_mesh(scalp, 'facecolor', 'skin')
lighting phong
camlight left
camlight right
material dull
alpha 0.3
ft_plot_mesh(headshape.fid, 'vertexcolor', 'k', 'vertexsize', 10);

%% Optional code: Check coregistration across files

figure
ft_plot_headmodel(headmodel_meg);
hold on
ft_plot_sens(grad, 'style', '*b');
ft_plot_headshape(headshape);

figure
ft_plot_headmodel(headmodel_meg);
hold on
ft_plot_sens(grad2, 'style', '*b');
ft_plot_headshape(headshape);

figure
ft_plot_headmodel(headmodel_meg);
hold on
ft_plot_sens(grad3, 'style', '*b');
ft_plot_headshape(headshape);

%% Clean-Up
rmpath(fullfile('..','..'))
