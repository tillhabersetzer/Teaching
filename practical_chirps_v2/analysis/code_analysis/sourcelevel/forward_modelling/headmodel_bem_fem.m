%--------------------------------------------------------------------------
% Within this script a BEM und FEM headmodels are computed.
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

%% BEM headmodel
%--------------------------------------------------------------------------
% https://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/

% Load coregistered mri
%----------------------
mri_segmented = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_T1w-segmented.mat']));

cfg              = [];
cfg.funparameter = 'scalp';
ft_sourceplot(cfg,mri_segmented);

% Mesh
%-----
cfg             = [];
cfg.tissue      = {'brain','skull','scalp'};
cfg.numvertices = [3000 2000 1000];
bnd             = ft_prepare_mesh(cfg,mri_segmented);

% Check surfaces
%---------------
figure;
ft_plot_mesh(bnd(1), 'facecolor','none','edgecolor', 'k'); 
figure;
ft_plot_mesh(bnd(2), 'facecolor','none','edgecolor', 'k'); 
figure;
ft_plot_mesh(bnd(3), 'facecolor','none','edgecolor', 'k'); 

% In case segmentation and surface meshes failed (holes, spikes etc.?!)
%--------------------------------------------------------------------------
% - bias correction /  Inhomogeneous anatomical image
%   https://www.fieldtriptoolbox.org/faq/why_does_my_eegheadmodel_look_funny/
% - if there are still spikes/ holes in meshes, change number of vertices

failure = 0;
if failure
    % Load_resliced mri
    %------------------
    mri_resliced  = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_T1w-resliced.mat']));

    % Bias correction
    %----------------
    mri_corrected = ft_volumebiascorrect([], mri_resliced);

    cfg              = [];
    ft_sourceplot(cfg,mri_corrected);
    
    % Segmentation
    %-------------
    cfg               = [];
    cfg.output        = {'brain','skull','scalp'}; 
    mri_corrected_seg = ft_volumesegment(cfg, mri_corrected);
    
    cfg = [];
    cfg.funparameter = 'scalp';
    ft_sourceplot(cfg,mri_corrected_seg);

    % Mesh
    %-----
    cfg             = [];
    cfg.method      = 'projectmesh';
    cfg.tissue      = {'brain','skull','scalp'};
    cfg.numvertices = [3000 2000 1000]; % change numbers
    bnd             = ft_prepare_mesh(cfg,mri_corrected_seg);

    % Other approach
    %---------------
    % https://mailman.science.ru.nl/pipermail/fieldtrip/2020-December/040547.html
    % https://github.com/meeg-cfin/nemolab/blob/master/basics/nemo_mriproc.m
    addpath(settings.path2iso2mesh) % iso2mesh
    cfg             = [];
    cfg.method      = 'projectmesh';
    cfg.tissue      = {'brain','skull','scalp'};
    cfg.numvertices = 10000;
%     cfg.numvertices = [10000,10000,1000];    % We'll decimate later
    cfg.spmversion = 'spm12';
    bnd            = ft_prepare_mesh(cfg, mri_corrected_seg);
    % Mesh repairs - Not yet implemented in FT
    %    targetsize = [500 1000 1500 1500]; % final mesh size desired for layers (in order given above)
    targetsize = [1000 1000 1000]; % final mesh size desired for layers (in order given above)
    
    % decimate, check, and repair individual meshes using iso2mesh
    for ii = 1:length(bnd)
        [bnd(ii).pos, bnd(ii).tri] = meshresample(bnd(ii).pos, bnd(ii).tri, targetsize(ii)/size(bnd(ii).pos,1));
        [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'dup');
        [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'isolated');
        [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'deep');
        [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'meshfix');
    end  
    %% Ensure no overlaps - seems not to help, dont use it
%     bnd = nemo_decouplesurf(bnd); % https://github.com/meeg-cfin/nemolab/blob/master/basics/nemo_decouplesurf.m

    % Check surfaces
    %---------------
    figure;
    ft_plot_mesh(bnd(1), 'facecolor','none','edgecolor', 'r'); % brain
    figure;
    ft_plot_mesh(bnd(2), 'facecolor','none','edgecolor', 'k'); % skull
    figure;
    ft_plot_mesh(bnd(3), 'facecolor','none','edgecolor', 'k'); % scalp

    % all surfaces together
    figure
    ft_plot_mesh(bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
    hold on;
    ft_plot_mesh(bnd(2),'edgecolor','none','facealpha',0.4);
    hold on;
    ft_plot_mesh(bnd(1),'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha',0.5); 

end

% Headmodel
%----------
% need to install Openmeeg in case method 'openmeeg' is used
setenv('PATH', fullfile(settings.path2openmeeg,'bin'))
setenv('LD_LIBRARY_PATH', fullfile(settings.path2openmeeg,'lib'))
cfg           = [];
cfg.method    = 'openmeeg'; % OpenMEEG needs to be installed
headmodel_bem = ft_prepare_headmodel(cfg, bnd);

cfg             = [];
cfg.method      = 'bemcp';
headmodel_bemcp = ft_prepare_headmodel(cfg, bnd);

% Visualization
%--------------
ft_plot_mesh(headmodel_bem.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(headmodel_bem.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(headmodel_bem.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);

ft_plot_mesh(headmodel_bemcp.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(headmodel_bemcp.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(headmodel_bemcp.bnd(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);

% Save headmodel
%---------------
save(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_headmodel-bem.mat']),'headmodel_bem'); % in mm
save(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_headmodel-bemcp.mat']),'headmodel_bemcp'); % in mm

%% FEM headmodel
%--------------------------------------------------------------------------
% https://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_fem/

% Load coregistered mri
%----------------------
% mri_resliced = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_T1w-resliced.mat']));
mri_coreg = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_T1w-coregistered.mat']));

% If needed: Bias correction
%---------------------------
% mri_resliced = ft_volumebiascorrect([], mri_resliced);
mri_coreg = ft_volumebiascorrect([], mri_coreg);

tissue_order          = {'gray','white','csf','skull','scalp'};
tissue_conductivities = [0.33 0.14 1.79 0.01 0.43]; % fieldtrip

cfg           = [];
cfg.output    = tissue_order; % conductivities: [0.33 0.14 1.79 0.01 0.43]
mri_segmented = ft_volumesegment(cfg, mri_coreg);
% mri_segmented = ft_volumesegment(cfg, mri_resliced);

% Check segmentation
%-------------------
seg_i = ft_datatype_segmentation(mri_segmented,'segmentationstyle','indexed');

cfg              = [];
cfg.funparameter = 'tissue';
cfg.funcolormap  = lines(6); % distinct color per tissue
cfg.location     = 'center';
cfg.atlas        = seg_i;    % the segmentation can also be used as atlas
ft_sourceplot(cfg, seg_i);

% Mesh
%-----
cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
mesh       = ft_prepare_mesh(cfg,mri_segmented);

% Visualize mesh
%---------------
addpath(fullfile('..','..','helper_functions'))
% csf: 1, gm: 2, scalp: 3, skull: 4, wm: 5
ts = 3; % tissuetype for plotting
figure
mesh2 =[];
mesh2.hex    = mesh.hex(mesh.tissue==ts,:); %mesh2.hex(1:size(mesh2.hex),:);
mesh2.pos    = mesh.pos;
mesh2.tissue = mesh.tissue(mesh.tissue==ts,:); %mesh.tissue(1:size(mesh2.hex),:);

mesh_ed = th_mesh2edge(mesh2);
patch('Faces',mesh_ed.poly,...
  'Vertices',mesh_ed.pos,...
  'FaceAlpha',.5,...
  'LineStyle','none',...
  'FaceColor',[1 1 1],...
  'FaceLighting','gouraud');

xlabel('coronal');
ylabel('sagital');
zlabel('axial')
camlight;
axis on;

% Visualization
%--------------
ft_plot_mesh(mesh, 'surfaceonly', 'yes');

% Headmodel
%----------
% Adapt order of conductivities
conductivities = zeros(1,5);
for t=1:5
    idx               = contains(tissue_order,mesh.tissuelabel{t});
    conductivities(t) = tissue_conductivities(idx);
end

cfg              = [];
% cfg.conductivity = [0.43 0.0024 1.79 0.14 0.33]; % same as tissuelabel in vol_simbio
% cfg.tissuelabel  = {'scalp', 'skull', 'csf', 'gray','white'};
cfg.method       ='simbio';
cfg.conductivity = conductivities;  % order follows mesh.tissuelabel
headmodel_fem    = ft_prepare_headmodel(cfg, mesh);

% Save headmodel and mesh
%------------------------
save(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_mesh-fem.mat']),'mesh'); % in mm
save(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_headmodel-fem.mat']),'headmodel_fem','-v7.3'); % in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bnd = nemo_decouplesurf(bnd)
for ii = 1:length(bnd)-1
    disp(num2str(ii))
  % Despite what the instructions for surfboolean says, surfaces should
  % be ordered from inside-out!! (ii=1: brain, ii=2: skull, ii=3: scalp)
  % smaller indices first: inside-out
  [newnode, newelem] = surfboolean(bnd(ii).pos,bnd(ii).tri,'decouple',bnd(ii+1).pos,bnd(ii+1).tri);
%   newelem = newelem(newelem(:,4)==2,1:3);
%   newnode = newnode(newnode(:,4)==2,1:3);
  
%  bnd(ii+1).tri = newelem(newelem(:,4)==2,1:3) - size(bnd(ii+1).pos,1);
%  bnd(ii+1).pos = newnode(newnode(:,4)==2,1:3);
  bnd(ii).tri = newelem; % - size(bnd(ii+1).pos,1);
  bnd(ii).pos = newnode;

%   plotmesh(bnd(3).pos,bnd(3).tri)
%   plotmesh(newnode,newelem)
  % cfg.tissue = {'brain','skull','scalp'};
end
end
