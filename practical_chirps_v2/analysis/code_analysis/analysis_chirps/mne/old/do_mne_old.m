close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!SCRIPT IS DESIGNED FOR LINUX!!!!
% to work on windows replace new = '/media/till/Samsung_T5' with new = old
% then everthing with the pays stays the same
%
% This script tries to give you a first overview about the averages and
% minimum norm estimates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject03'};

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [4, 30];

% apply ica
ica_status = 0;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% apply baseline correction
baseline_correction_status = 1;

% downsample data
downsample_status = 0;
fs_down           = 300;
%--------------------------------------------------------------------------

% addpath for subject_files information
old = 'Z:'; new = '/media/till/Samsung_T5';
new = 'Z:';

% addpath for subject_files information
addpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new));
% addpath for ica functions
addpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new));
% addpath for preprocessing function
addpath(replace(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions'],old,new));

%% analysis

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj

%% preprocessing
config.subject                    = subjects{s};
config.files2preproc              = files2preproc;
config.bp_freq                    = bp_freq;
config.ica_status                 = ica_status;
config.ica_files                  = ica_files;
config.baseline_correction_status = baseline_correction_status;
config.downsample_status          = downsample_status;
config.fs_down                    = fs_down;

data_preprocessed = preprocess_chirps(config);

eval(subjects{s})
% change paths if you work with Linux
subjectdata.rawdatadir = replace(subjectdata.rawdatadir,old,new);    
    
%% Averaging and noise-covariance estimation
% for a correct noise-covariance estimation it is important that 
% you used the cfg.demean = 'yes';

coi                  = 'meg';

cfg                  = [];
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = [-inf 0]; % timepoints before the zero timepoints in trials
avg                  = ft_timelockanalysis(cfg, data_preprocessed);

% estimate noise covariance different (empty rooms)
empty_room = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
N          = length(empty_room);
noise      = cell(1,N);

for n = 1:N
    cfg              = [];
    cfg.dataset      = replace([subjectdata.rawdatadir filesep empty_room{n} '.fif'],old,new);
    cfg.channel      = coi; 
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;            % ensure that sensors are expressed in SI units
    noise{n}         = ft_preprocessing(cfg);   
end

hdr                = noise{1}.hdr;
cfg                = [];
cfg.keepsampleinfo = 'no';
noise              = ft_appenddata(cfg,noise{:});
noise.hdr          = hdr;

cfg                  = [];
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg,noise);

figure
subplot(1,2,1)
imagesc(avg_noise.cov);
title('empty room measurements')
subplot(1,2,2)
imagesc(avg.cov);
title('baseline interval')
  
%% Forward solution
% load headmodel
headmodel = importdata(replace(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_headmodel.mat']),old,new)); % mm
% load sourcemodel
sourcemodel = importdata(replace(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_4k.mat']),old,new)); % mm

% prepare leadfield
%------------------
cfg                       = [];
cfg.grad                  = data_preprocessed.grad; % sensor information
cfg.channel               = coi;                    % the used channels
cfg.sourcemodel           = sourcemodel;            % source points
cfg.headmodel             = headmodel;              % volume conduction model
cfg.singleshell.batchsize = 5000;                   % speeds up the computation
leadfield                 = ft_prepare_leadfield(cfg);

% check
%------
% same units for plotting
dataset = replace([subjectdata.rawdatadir filesep 'fixer_olsa_tsss_mc.fif'],old,new);
shape   = ft_read_headshape(dataset,'unit','mm');
grad    = ft_convert_units(ft_read_sens(dataset,'senstype','meg'),'mm'); 

figure
ft_plot_headmodel(headmodel,'facealpha',0.5);
hold on
ft_plot_sens(grad, 'style', '*b');
ft_plot_headshape(shape);
ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', ...
             'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);
  
%% Inverse solution
% modify noise covariance estimation
avg.cov = avg_noise.cov;

cfg                    = [];
cfg.method             = 'mne';
cfg.sourcemodel        = leadfield;
cfg.headmodel          = headmodel;
cfg.mne.prewhiten      = 'yes'; % combine different sensortypes with different units
cfg.mne.lambda         = 3; % scaling factor for noise
cfg.mne.scalesourcecov = 'yes';
cfg.mne.keepfilter     = 'yes';
source                 = ft_sourceanalysis(cfg,avg);

%source.inside = mask;
source.coordsys = 'neuromag';
%% Visualize
signal = source.avg.pow(:,640); % 100ms

figure
ft_plot_mesh(source, 'vertexcolor', signal);
view([180 0]);
h = light; 
set(h, 'position', [0 1 0.2]); 
lighting gouraud;
material dull;

% movie
cfg            = [];
%cfg.projectmom = 'yes';
sd             = ft_sourcedescriptives(cfg,source);

cfg              = [];
cfg.atlas        = atlas;
cfg.roi          = {'R_superiortemporal'};
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,sd);

cfg              = [];
%cfg.maskparameter = mask;
cfg.atlas        = atlas;
cfg.roi          = {'R_superiortemporal'};
cfg.funparameter = 'pow';
cfg.method       = 'surface';
cfg.latency      = 0.140;
ft_sourceplot(cfg,source)

surf         = source;
surf.avg.pow = source.avg.pow(:,640);
cfg          = [];
handles = myownsurfplot(cfg,sourcemodel,surf);

%% Inverse solution with atlas
atlas_orig = ft_read_atlas('Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\template\atlas\aal\ROI_MNI_V4.nii');
%path = replace(fullfile(subjectdata.sourcemodel,'workbench',[subjectdata.subjectname,'.L.aparc.a2009s.4k_fs_LR.label.gii']),old,new);
path_label_r = replace(fullfile(subjectdata.sourcemodel,'workbench',[subjectdata.subjectname,'.R.aparc.4k_fs_LR.label.gii']),old,new);
path_mesh_r  = replace(fullfile(subjectdata.sourcemodel,'workbench',[subjectdata.subjectname,'.R.midthickness.4k_fs_LR.surf.gii']),old,new);

atlas_r = ft_read_atlas({path_label_r,path_mesh_r});
atlas_l = ft_read_atlas({strrep(path_label_r,'.R.','.L.'),strrep(path_mesh_r,'.R.','.L.')});

atlas.parcellation = [atlas_l.parcellation;atlas_r.parcellation];
atlas.parcellationlabel = [atlas_l.parcellationlabel;atlas_r.parcellationlabel];
atlas.rgba = [atlas_l.rgba;atlas_r.rgba];
atlas.pos = [atlas_l.pos;atlas_r.pos];
atlas.tri = [atlas_l.tri;atlas_r.tri];
atlas.unit = 'mm';
atlas.coordsys = 'neuromag';

cfg       = [];
cfg.roi   = {'R_superiortemporal'};
cfg.atlas = atlas;
cfg.inputcoord = 'unknown';
mask      = ft_volumelookup(cfg,sourcemodel);




atlas_l = ft_read_atlas(path); 
atlas_r = ft_read_atlas(strrep(path,'.L.','.R.'));

roi                     = {'L_superiortemporal','R_superiortemporal'};
mask                    = mne_give_maks(roi,atlas_l,atlas_r,sourcemodel);
sourcemodel_mask        = sourcemodel;
sourcemodel_mask.inside = mask;

% prepare leadfield
%------------------
cfg                       = [];
cfg.grad                  = data_preprocessed.grad; % sensor information
cfg.channel               = coi;                    % the used channels
cfg.sourcemodel           = sourcemodel_mask;       % source points
cfg.headmodel             = headmodel;              % volume conduction model
cfg.singleshell.batchsize = 5000;                   % speeds up the computation
leadfield_mask            = ft_prepare_leadfield(cfg);


% modify noise covariance estimation
avg.cov = avg_noise.cov;

cfg                    = [];
cfg.method             = 'mne';
cfg.sourcemodel        = leadfield_mask;
cfg.headmodel          = headmodel;
cfg.mne.prewhiten      = 'yes'; % combine different sensortypes with different units
cfg.mne.lambda         = 3; % scaling factor for noise
cfg.mne.scalesourcecov = 'yes';
cfg.mne.keepfilter     = 'yes';
source_mask            = ft_sourceanalysis(cfg,avg);
%% Visualize with atlas
signal = source_mask.avg.pow(:,640); % 100ms

figure
ft_plot_mesh(source_mask, 'vertexcolor', signal);
view([180 0]);
h = light; 
set(h, 'position', [0 1 0.2]); 
lighting gouraud;
material dull;

% movie
cfg              = [];
%cfg.projectmom = 'yes';
sd_mask          = ft_sourcedescriptives(cfg,source_mask);

cfg              = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,sd_mask);

end

%% Clean up
rmpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new))
rmpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new))
rmpath(replace(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions'],old,new));

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,avg);
suptitle('average magnetometer')

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
%cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,avg);
suptitle('average gradiometer')

% atlas information
%roi = ft_read_headshape('/media/till/Volume/mne_workbench/template_files/R.atlasroi.4k_fs_LR.shape.gii');
%surf = ft_read_headshape('/media/till/Volume/mne_workbench/template_files/R.sphere.4k_fs_LR.surf.gii');

%label = ft_read_headshape('/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.L.aparc.4k_fs_LR.label.gii');
%label1 = ft_read_atlas('/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.R.aparc.4k_fs_LR.label.gii');
%label2 = ft_read_atlas('/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.R.aparc.a2009s.4k_fs_LR.label.gii');
%label = ft_read_atlas('/media/till/Volume/mne_workbench/template_files/R.atlasroi.4k_fs_LR.shape.gii');

% a = '/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.R.aparc.a2009s.4k_fs_LR.label.gii';
% b = '/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.R.aparc.4k_fs_LR.label.gii';
% label_a = ft_read_atlas({a, strrep(a,'.R.','.L.')});
% label_b = ft_read_atlas({b, strrep(b,'.R.','.L.')});


function [mask] = mne_give_maks(roi,atlas_l,atlas_r,sourcemodel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a brain mask for a given sourcemodel based on atlases for each
% hemisphere.
% - hcp workbench atlases
% - mne sourcemodel
% - region of interests (roi) defined as a string or cell array of strings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do it for each hemisphere separately

% find correct parcellation number
label_atlas_l = find(contains(atlas_r.parcellationlabel,roi,'IgnoreCase',true));
label_atlas_r = find(contains(atlas_l.parcellationlabel,roi,'IgnoreCase',true));

% find correct surface points for given parcellation number
idx_atlas_l = ismember(atlas_l.parcellation,label_atlas_l);
idx_atlas_r = ismember(atlas_r.parcellation,label_atlas_r);

% check for left/right hemisphere order
label_hemisphere_l = find(contains(sourcemodel.brainstructurelabel,'Left','IgnoreCase',true));
label_hemisphere_r = find(contains(sourcemodel.brainstructurelabel,'Right','IgnoreCase',true));

% find correct hemisphere indices in sourcemodel
idx_hemisphere_l = find(ismember(sourcemodel.brainstructure,label_hemisphere_l));
idx_hemisphere_r = find(ismember(sourcemodel.brainstructure,label_hemisphere_r));

% create mask
mask = zeros(size(sourcemodel.inside));
mask(idx_hemisphere_l) = idx_atlas_l;
mask(idx_hemisphere_r) = idx_atlas_r;
mask = logical(mask);
end

function [handles] = myownsurfplot(cfg,surf,source)

% function [handles] = myownsurfplot(cfg,surf,source)
%
% allows to plot source results obtained from a cortical mesh
%
% cfg - opacity: specifies the opacities of the vertices e.g. [0,1]
% cfg - colormapping: defines the colorrange of the vertices e.g. [0,1]
% cfg - edgealpha: defines the transparency of the edges for identifying gyri & sulci, e.g. 0.2
% cfg - mask matrix specifying opacity, format: [n,1] for n vertices

if nargin < 3
  cfg = struct();
end


cortex_light = [0.781 0.762 0.664];
cortex_dark  = [0.781 0.762 0.664]/2;


sourcevals = source.avg.pow(:);
backgcolor = repmat(cortex_light, size(surf.pos,1), 1);

if ~isfield(cfg,'opacity')
  opacmin = min((source.avg.pow(:)));
  opacmax = max((source.avg.pow(:)));
else
  opacmin = cfg.opacity(1);
  opacmax = cfg.opacity(2);
end
if ~isfield(cfg,'colormapping')
  fcolmin = min((source.avg.pow(:)));
  fcolmax = max((source.avg.pow(:)));
else
  fcolmin = cfg.opacity(1);
  fcolmax = cfg.opacity(2);
end
if ~isfield(cfg,'edgealpha'), edgealpha = 1; end
if ~isfield(cfg,'mask')
  maskval = source.avg.pow(:);
else
  maskval = cfg.mask;
end


figure;

h1 = patch('Vertices', surf.pos, 'Faces', surf.tri, 'FaceVertexCData', backgcolor , 'FaceColor', 'interp');
%set(h1, 'EdgeColor', 'none');
set(h1, 'EdgeColor', [0,0,0],'EdgeAlpha',edgealpha);
axis   off;
axis vis3d;
axis equal;

h2 = patch('Vertices', surf.pos, 'Faces', surf.tri, 'FaceVertexCData', sourcevals , 'FaceColor', 'interp');
%set(h2, 'EdgeColor',  'none');
set(h2, 'EdgeColor', [0,0,0],'EdgeAlpha',edgealpha);
set(h2, 'FaceVertexAlphaData', maskval);
set(h2, 'FaceAlpha',          'interp');
set(h2, 'AlphaDataMapping',   'scaled');
alim(gca, [opacmin opacmax]);
caxis(gca,[fcolmin fcolmax]);
lighting gouraud

colorbar;

handles.h1 = h1;
handles.h2 = h2;

end
