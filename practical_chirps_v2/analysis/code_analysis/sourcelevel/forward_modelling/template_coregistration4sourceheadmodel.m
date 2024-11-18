%--------------------------------------------------------------------------
% Till Habersetzer (21.01.2022)
% This script is designed to handle template mri's in case a participant
% has no recorded mri. 
% Therefore, this script is quite extensive and serves multiple purposes. 
% 1.) Coregistration of template mri (here: MNI-ICBM152) and recorded
%     headshape during digitization
% 2.) Creation of a singleshell headmodel based on coregistered mri
% 3.) Creation of a 3D-grid-based volume conduction model for beamforming
% 
% Step 2.) and 3.) are based on the following script and contain same code.
% - headmodel.m
% - volumetric_sourcemodel.m
% Step 1.) comprises two successively applied coordinate transformations.
% The first alignement (t1) is rather coarse and the second alignment (t2)
% is a subsequent refinement. 
% The idea is based on: https://www.fieldtriptoolbox.org/example/fittemplate/
% In this script only a coregistered headmodel is created but no mri or 
% sourcemodel. The coregistration is based on a transformation between the
% digitized headshape and an already precomputed headmodel in fieldtrip
% which is subsequently transformed according to the digitized headshape.
% Here, in contrast, the coregistration is based on a template mri from
% which the headshape is extracted and coregistered to the digitized
% headshape. Hence, the template mri itself can be coregistred to the 
% neuromag coordinate system based on the subjects digitized headshape.
% Following this, headmodel and sourcemodel can be computed from the
% coregistred mri. 
%
% Remark: The initial idea for the creation of a coregistered headmodel
% based on the fieldtrip webpage is implemented in: template_headmodel.m 
% -------------------------------------------------------------------------

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
subidx  = 1;
subject = ['sub-',num2str(subidx,'%02d')];


% choose dataset for headshape extraction
megfile  = fullfile(settings.path2bids,subject,'meg',[subject,'_task-clicks_meg.fif']); 
dir2save = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling');

% addpath for helper_functions
%-----------------------------
addpath(fullfile('..','..','helper_functions'))
%--------------------------------------------------------------------------

%% Coregistration between template MRI and Neuromag
% Determine coordinate transformation between template mri (MNI-ICBM152) and neuromag headhshape
%------------------------------------------------------------------------------------------------

% Load Data
%--------------------------------------------------------------------------

% Extract headshape
headshape = ft_read_headshape(megfile,'unit','mm');

mri          = ft_read_mri(fullfile(settings.path2bids,'derivatives','mni_icbm152_nlin_sym_09a','mni_icbm152_t1_tal_nlin_sym_09a.nii')); % mm
mri.coordsys = 'mni';

% Extract headshape from template
%--------------------------------------------------------------------------
cfg                = [];
cfg.output         = {'scalp'}; 
cfg.scalpthreshold = 0.2; % Remove Aliasing artifact at the top that leads to horns, see https://www.fieldtriptoolbox.org/faq/why_does_my_eegheadmodel_look_funny/
cfg.scalpsmooth    = 'no';
mri_segmented      = ft_volumesegment(cfg, mri);

mri_segmented.anatomy = mri.anatomy;

cfg = [];
ft_sourceplot(cfg, mri_segmented);

% Create mesh
cfg             = [];
cfg.tissue      = 'scalp';
cfg.numvertices = 2000;
template        = ft_prepare_mesh(cfg,mri_segmented);

figure;
ft_plot_mesh(template,'facecolor','none'); % scalp

% Coregistration
%--------------------------------------------------------------------------
% coregister the template headshape with the Polhemus head shape

cfg                      = [];
cfg.template.headshape   = headshape;
cfg.checksize            = inf;
cfg.individual.headshape = template;
cfg                      = ft_interactiverealign(cfg);

% save transformation matrix and apply to template headshape and template mri
t1          = cfg.m;
template_t1 = ft_transform_geometry(cfg.m, template);
mri_t1      = ft_transform_geometry(cfg.m, mri);

% Refining the transformation
%--------------------------------------------------------------------------

% Method 2: On the basis of the full head surface
%------------------------------------------------
% This requires an external toolbox
addpath(genpath(settings.path2iso2cpd2)) 

% It is important that the template scalp surface only contains features 
% that are also in the Polhemus surface, and vice versa. We can use 
% ft_defacemesh to remove some features.

% Visualize both meshes:
figure;
ft_plot_mesh(template_t1);  
ft_plot_mesh(headshape); 

% Remove details so that template and headshape are similar in terms of
% head converage
defaced_template = template_t1;

cfg              = [];
cfg.translate    = [-40 0 -50];
cfg.scale        = [200 200 200];
cfg.rotate       = [0 0 0];
defaced_template = ft_defacemesh(cfg, defaced_template);

cfg               = [];
cfg.translate     = [-40 0 -50];
cfg.scale         = [200 200 200];
cfg.rotate        = [0 0 0];
defaced_headshape = ft_defacemesh(cfg, headshape);

% Have another look how well the surfaces match
figure;
ft_plot_mesh(defaced_template);  
ft_plot_mesh(defaced_headshape); 

% Determine the transformation and apply it to all 3 surfaces of the template head model
% cfg             = [];
% cfg.headshape   = defaced_headshape;
% cfg.template    = defaced_template;
% cfg.method      = 'fittemplate';
% [template_fit_t2,cfg] = ft_prepare_mesh(cfg, template_t1);

% -> replace with lower level function to retrieve transformation matrix 
% see ft_prepare_mesh (line 200-205)
% prepare_mesh_fittemplate is in private folder and cant be accessed 
% through addpath, therefore its copied
t2          = prepare_mesh_fittemplate(defaced_headshape.pos, defaced_template.pos);
template_t2 = ft_transform_geometry(t2, template_t1);
mri_t2      = ft_transform_geometry(t2, mri_t1);

% Have another look how well the surfaces match
figure;
ft_plot_mesh(template_t2);  
ft_plot_mesh(headshape); 

cfg = [];
ft_sourceplot(cfg, mri_t2);

% Save data
%----------
if ~exist(fullfile(dir2save,'headmodel'), 'dir')
    mkdir(fullfile(dir2save,'headmodel'))
end
save(fullfile(dir2save,'headmodel',[subject,'_template-coregistration','.mat']),'t1','t2','mri_t1','mri_t2','template_t1','template_t2'); 

%% Compute headmodel based on template mri
%--------------------------------------------------------------------------

% Depending on the coregistration mri_t1 and mri_t2 can be used as a
% startung point. Even though, mri_t2 is a refinement of mri_t1, it might
% be less suited because of the automatic alignment algorithm which may
% lead to suspicious and faulty results. Check the sourceplots after the
% reslicing!
mri_coreg = importdata(fullfile(dir2save,'headmodel',[subject,'_template-coregistration','.mat']));
mri_t2    = mri_coreg.mri_t2;
mri_t1    = mri_coreg.mri_t1;
clear mri_coreg;

% 0. Coregistration
%------------------
% use coregistered mri
% mri_orig = mri_t2;
mri_orig = mri_t1;

% 1. Reslice
%-----------
cfg            = [];
cfg.resolution = 1;
cfg.dim        = [256 256 256];
cfg.zrange     = [-100,155]; % could be adjusted, see above remark, sum up to 255 -> 256 voxel
cfg.coordsys   = 'neuromag';
mri_resliced   = ft_volumereslice(cfg,mri_orig);

cfg = [];
ft_sourceplot(cfg, mri_orig);
cfg = [];
ft_sourceplot(cfg, mri_resliced);

% Save the Resliced Anatomy in a FreeSurfer Compatible Format
%------------------------------------------------------------
% Volumewrite saves no Information about the Coordinate System!
cfg           = [];
cfg.filename  = fullfile(dir2save,'headmodel',[subject,'_T1w-coregistered.mgz']);
cfg.filetype  = 'mgz';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_resliced);
% --> Needed for Headmodel and FreeSurfer

% 2. Segmentation (most time consuming)
%--------------------------------------
cfg           = [];
cfg.output    = {'brain','skull','scalp'}; % you only need brain for meg
mri_segmented = ft_volumesegment(cfg, mri_resliced);

% Copy the Anatomy into the segmented Mri
mri_segmented.anatomy = mri_resliced.anatomy;

% Check
%------
cfg              = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, mri_segmented);

% 3. Headmodel
%-------------
cfg        = [];
cfg.method = 'singleshell';
vol        = ft_prepare_headmodel(cfg, mri_segmented);

save(fullfile(dir2save,'headmodel',[subject,'_T1w-resliced.mat']),'mri_resliced'); % in mm
save(fullfile(dir2save,'headmodel',[subject,'_T1w-segmented.mat']),'mri_segmented'); % in mm
save(fullfile(dir2save,'headmodel',[subject,'_headmodel-singleshell.mat']),'vol'); % in mm

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

%% Compute beamformer sourcemodel based on template mri
%--------------------------------------------------------------------------

% load template grid
%-------------------
template_grid = importdata(fullfile(settings.path2fieldtrip,'template','sourcemodel','standard_sourcemodel3d10mm.mat'));
template_grid = ft_convert_units(template_grid,'mm'); 

% plot the atlas based grid
figure
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

% Coregistration of subject specific grid to the atlas based template grid
%-------------------------------------------------------------------------

% make individual subject's grid
mri_segmented = importdata(fullfile(dir2save,'headmodel',[subject,'_T1w-segmented.mat'])); % mm
% mri_segmented = ft_convert_units(mri_segmented,'m'); % use SI units

cfg           = [];
cfg.warpmni   = 'yes';
cfg.template  = template_grid;
cfg.nonlinear = 'yes';
cfg.mri       = mri_segmented;
cfg.unit      = 'mm';
sourcemodel   = ft_prepare_sourcemodel(cfg);

% Save data
%----------
if ~exist(fullfile(dir2save,'sourcemodel'), 'dir')
    mkdir(fullfile(dir2save,'sourcemodel'))
end

save(fullfile(dir2save,'sourcemodel',[subject,'_sourcemodel-volumetric.mat']),'sourcemodel'); % in mm

%% Cortical Sheet based sourcemodel
%--------------------------------------------------------------------------

check = 1;

% resolution - specify only one hemisphere
resolution{1} = '.L.midthickness.4k_fs_LR.surf.gii'; label{1} = '4k';
resolution{2} = '.L.midthickness.8k_fs_LR.surf.gii'; label{2} = '8k';

for r = 1:length(resolution)
        
    filename = [subject,resolution{r}];
    path     = fullfile(settings.path2bids,'derivatives',subject,'freesurfer','workbench',filename);
    
    % add right hemisphere
    files       = {path, strrep(path,'.L.','.R.')};
    sourcemodel = ft_read_headshape(files); 
    
    sourcemodel          = ft_determine_units(sourcemodel);
    sourcemodel.coordsys = 'neuromag';
    % sourcemodel.inside   = sourcemodel.atlasroi>0; -> seems not to work with ft_sourcestatistics !
    % sourcemodel          = rmfield(sourcemodel,'atlasroi');

    % add inflated surface model for visualization
    inflated          = ft_read_headshape(strrep(files, 'midthickness', 'inflated'));
    inflated          = ft_determine_units(inflated);
    inflated.coordsys = 'neuromag';
    % inflated.inside   = inflated.atlasroi>0;
    % inflated          = rmfield(inflated,'atlasroi');

    if check
        % check sourcemodel 
        megfile   = fullfile(settings.path2bids,subject,'meg',[subject,'_task-clicks_meg.fif']);
        shape     = ft_read_headshape(megfile,'unit','mm');
        grad      = ft_convert_units(ft_read_sens(megfile,'senstype','meg'),'mm'); 
        headmodel = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_headmodel-singleshell.mat'])); % mm
     
        figure
        ft_plot_headmodel(headmodel,'facealpha',0.1);
        hold on
        ft_plot_sens(grad, 'style', '*b');
        ft_plot_headshape(shape);
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', ...
                     'facealpha', 0.25, 'edgecolor', 'red','edgeopacity', 0.3); 
        title([subject,': ',label{r}])
    end

   % Save data
   %----------
   if ~exist(fullfile(dir2save,'sourcemodel'), 'dir')
    mkdir(fullfile(dir2save,'sourcemodel'))
   end

    save(fullfile(dir2save,'sourcemodel',[subject,'_sourcemodel-corticalsheet',label{r},'.mat']),'sourcemodel'); % in mm
    save(fullfile(dir2save,'sourcemodel',[subject,'_sourcemodel-inflatedcorticalsheet',label{r},'.mat']),'inflated'); % in mm

end % resolutions





