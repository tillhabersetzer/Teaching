%--------------------------------------------------------------------------
% checkout https://www.fieldtriptoolbox.org/workshop/baci2017/forwardproblem/
% https://www.fieldtriptoolbox.org/example/fem/
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
megfile = fullfile(settings.path2bids,subject,'meg',[subject,'_task-clicks_meg.fif']); 
%--------------------------------------------------------------------------

%% Load data
%--------------------------------------------------------------------------
elec    = ft_read_sens(megfile,'senstype','eeg'); % cm
grad    = ft_read_sens(megfile,'senstype','meg'); % cm
elec    = ft_convert_units(elec,'m');
grad    = ft_convert_units(grad,'m');

path2headmodel = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel');
% headmodel_bem = importdata(fullfile(path2headmodel,[subject,'_headmodel-bemcp.mat'])); % mm bemcp
headmodel_bem  = importdata(fullfile(path2headmodel,[subject,'_headmodel-bem.mat'])); % mm openmeeg
headmodel_fem  = importdata(fullfile(path2headmodel,[subject,'_headmodel-fem.mat'])); % mm
headmodel_meg  = importdata(fullfile(path2headmodel,[subject,'_headmodel-singleshell.mat'])); % mm

path2sourcemodel = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','sourcemodel');
sourcemodel_vol  = importdata(fullfile(path2sourcemodel,[subject,'_sourcemodel-volumetric.mat']));
sourcemodel_surf = importdata(fullfile(path2sourcemodel,[subject,'_sourcemodel-corticalsheet4k.mat']));

dir2save = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','leadfield');


%% Convert data in SI units
%--------------------------------------------------------------------------
headmodel_bem = ft_convert_units(headmodel_bem,'m');
headmodel_fem = ft_convert_units(headmodel_fem,'m');
headmodel_meg = ft_convert_units(headmodel_meg,'m');
          
sourcemodel_vol  = ft_convert_units(sourcemodel_vol,'m');
sourcemodel_surf = ft_convert_units(sourcemodel_surf,'m');

%% Computation of leadfields
%--------------------------------------------------------------------------

% MEG
%--------------------------------------------------------------------------
cfg                       = [];
cfg.sourcemodel           = sourcemodel_vol;
cfg.headmodel             = headmodel_meg;
cfg.grad                  = grad;
cfg.channel               = 'meg';
cfg.singleshell.batchsize = 'all'; % number of dipoles for which the leadfield is computed in a single call
leadfield_meg_vol         = ft_prepare_leadfield(cfg);

cfg.grid                  = sourcemodel_surf;
leadfield_meg_surf        = ft_prepare_leadfield(cfg);

% EEG
%--------------------------------------------------------------------------

% BEM
%----
% need to install Openmeeg in case method 'openmeeg' is used
setenv('PATH', fullfile(settings.path2openmeeg,'bin'))
setenv('LD_LIBRARY_PATH', fullfile(settings.path2openmeeg,'lib'))
cfg                = [];
cfg.sourcemodel    = sourcemodel_vol;
cfg.headmodel      = headmodel_bem;
cfg.elec           = elec;
cfg.channel        = 'eeg';
leadfield_bem_vol  = ft_prepare_leadfield(cfg);

cfg.sourcemodel     = sourcemodel_surf;
leadfield_bem_surf = ft_prepare_leadfield(cfg);

% FEM
%----
cfg                = [];
cfg.sourcemodel    = sourcemodel_vol;
cfg.headmodel      = headmodel_fem;
cfg.elec           = elec;
cfg.channel        = 'eeg';
leadfield_fem_vol  = ft_prepare_leadfield(cfg);

cfg.sourcemodel    = sourcemodel_surf;
leadfield_fem_surf = ft_prepare_leadfield(cfg);

%% Save data
%--------------------------------------------------------------------------

% make folder for data
%---------------------
if ~exist(dir2save, 'dir')
    mkdir(dir2save)
end

save(fullfile(dir2save,[subject,'_leadfield-meg-vol.mat']),'leadfield_meg_vol'); 
save(fullfile(dir2save,[subject,'_leadfield-meg-surf.mat']),'leadfield_meg_surf'); 

save(fullfile(dir2save,[subject,'_leadfield-bem-vol.mat']),'leadfield_bem_vol'); 
save(fullfile(dir2save,[subject,'_leadfield-bem-surf.mat']),'leadfield_bem_surf'); 

% save(fullfile(dir2save,[subject,'_leadfield-bemcp-vol.mat']),'leadfield_bem_vol'); 
% save(fullfile(dir2save,[subject,'_leadfield-bemcp-surf.mat']),'leadfield_bem_surf'); 

save(fullfile(dir2save,[subject,'_leadfield-fem-vol.mat']),'leadfield_fem_vol'); 
save(fullfile(dir2save,[subject,'_leadfield-fem-surf.mat']),'leadfield_fem_surf'); 

