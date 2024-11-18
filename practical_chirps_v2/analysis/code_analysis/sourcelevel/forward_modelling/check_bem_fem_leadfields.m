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
shape   = ft_read_headshape(megfile,'unit','m');
elec    = ft_read_sens(megfile,'senstype','eeg'); % cm
grad    = ft_read_sens(megfile,'senstype','meg'); % cm
elec    = ft_convert_units(elec,'m');
grad    = ft_convert_units(grad,'m');

path2headmodel = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel');
headmodel_bem  = importdata(fullfile(path2headmodel,[subject,'_headmodel-bem.mat'])); % mm openmeeg
% headmodel_bemcp = importdata(fullfile(path2headmodel,[subject,'_headmodel-bemcp.mat'])); % mm bemcp
headmodel_fem = importdata(fullfile(path2headmodel,[subject,'_headmodel-fem.mat'])); % mm

path2sourcemodel = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','sourcemodel');
sourcemodel_vol  = importdata(fullfile(path2sourcemodel,[subject,'_sourcemodel-volumetric.mat']));
sourcemodel_surf = importdata(fullfile(path2sourcemodel,[subject,'_sourcemodel-corticalsheet4k.mat']));

path2leadfield = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','leadfield');
lf_meg_vol  = importdata(fullfile(path2leadfield,[subject,'_leadfield-meg-vol.mat']));
lf_meg_surf = importdata(fullfile(path2leadfield,[subject,'_leadfield-meg-surf.mat']));

lf_bem_vol  = importdata(fullfile(path2leadfield,[subject,'_leadfield-bem-vol.mat'])); % openmeeg
lf_bem_surf = importdata(fullfile(path2leadfield,[subject,'_leadfield-bem-surf.mat'])); % openmeeg

% lf_bemcp_vol  = importdata(fullfile(path2leadfield,[subject,'_leadfield-bemcp-vol.mat'])); % bemcp
% lf_bemcp_surf = importdata(fullfile(path2leadfield,[subject,'_leadfield-bemcp-surf.mat'])); % bemcp
 
lf_fem_vol  = importdata(fullfile(path2leadfield,[subject,'_leadfield-fem-vol.mat']));
lf_fem_surf = importdata(fullfile(path2leadfield,[subject,'_leadfield-fem-surf.mat']));


headmodel_bem = ft_convert_units(headmodel_bem,'m');
headmodel_fem = ft_convert_units(headmodel_fem,'m');

sourcemodel_vol  = ft_convert_units(sourcemodel_vol,'m');
sourcemodel_surf = ft_convert_units(sourcemodel_surf,'m');

%% Visualization geometries and headmodels
%--------------------------------------------------------------------------

% BEM
%----
figure
% headmodel
ft_plot_mesh(headmodel_bem.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(headmodel_bem.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(headmodel_bem.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha',0.5); 
% sourcemodel
hold on
ft_plot_mesh(sourcemodel_surf, 'maskstyle', 'opacity', 'facecolor', 'black', ...
                         'facealpha', 0.25, 'edgecolor', 'black','edgeopacity', 0.3); 
% ft_plot_mesh(sourcemodel_vol.pos(sourcemodel_vol.inside,:)); % plot only locations inside the volume
% sensors shape
ft_plot_sens(grad, 'style', '*b');
ft_plot_sens(elec, 'elecshape', 'disc');
ft_plot_headshape(shape);

% FEM
%----
addpath(fullfile('..','..','helper_functions'))
% csf: 1, gm: 2, scalp: 3, skull: 4, wm: 5
ts = 3; % tissuetype for plotting
figure
mesh        = [];
mesh.hex    = headmodel_fem.hex(headmodel_fem.tissue==ts,:); %mesh2.hex(1:size(mesh2.hex),:);
mesh.pos    = headmodel_fem.pos;
mesh.tissue = headmodel_fem.tissue(headmodel_fem.tissue==ts,:); %mesh.tissue(1:size(mesh2.hex),:);

mesh_ed = th_mesh2edge(mesh);
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
ft_plot_sens(elec, 'elecshape', 'disc','facecolor','g');
ft_plot_sens(grad, 'style', '*b');
% sourcemodel
hold on
ft_plot_mesh(sourcemodel_surf, 'maskstyle', 'opacity', 'facecolor', 'black', ...
                         'facealpha', 0.25, 'edgecolor', 'black','edgeopacity', 0.3); 

%% Visualization leadfields - Dipole geometry
%--------------------------------------------------------------------------

% sourcemodel_type = 'vol';
sourcemodel_type = 'surf';

switch sourcemodel_type
    case 'vol'
        roi         = 'Heschl_R'; % roi for dipole
        sourcemodel = sourcemodel_vol;

        % Apply atlas
        template_grid = importdata(fullfile(settings.path2fieldtrip,'template','sourcemodel','standard_sourcemodel3d10mm.mat'));
        template_grid = ft_convert_units(template_grid,'mm'); 
        template_grid.coordsys = 'mni';
        atlas = ft_read_atlas(fullfile(settings.path2fieldtrip,'template','atlas','aal','ROI_MNI_V4.nii')); % mm
        
        cfg            = [];
        cfg.atlas      = atlas;
        cfg.roi        = roi;
        cfg.inputcoord = 'mni';
        mask           = ft_volumelookup(cfg,template_grid);

        idx            = find(mask);
        idxvol         = idx(1); % take first point for leadfield
        idxp           = idxvol;
    case 'surf'
        roi         = 'R_G_temp_sup-G_T_transv'; % roi for dipole
        sourcemodel = sourcemodel_surf;

        % Choose atlas
        % atlasname  = '.L.aparc.4k_fs_LR.label.gii';
        atlasname  = '.L.aparc.a2009s.4k_fs_LR.label.gii';
        atlas      = generate_surfatlas(atlasname,sourcemodel,settings.path2bids,subject);

        % pick a random location of the auditory cortex
        idx1      = find(strcmp(atlas.parcellationlabel,roi));
        idx2      = find(atlas.parcellation == idx1);
        % take first point for leadfield
        idxsurf   = idx2(1);
        idxp      = idxsurf;
end

% Visualize geometry
%-------------------
figure('Position', [1 1 1080 1080])
% Plot BEM headmodel
ft_plot_mesh(headmodel_bem.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.2, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(headmodel_bem.bnd(2),'edgecolor','none','facealpha',0.4);
ft_plot_mesh(headmodel_bem.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha',0.3);
% Plot sourcemodel
switch sourcemodel_type
    case 'vol'
        ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
    case 'surf'
        ft_plot_mesh(sourcemodel,'facecolor',[0.2 0.2 0.2], 'facealpha', 0.2, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
end
% Plot sensors
ft_plot_sens(grad,'unit',sourcemodel.unit,'style', '*b');
ft_plot_sens(elec,'unit',sourcemodel.unit);
% Add dipoles
ft_plot_dipole(sourcemodel.pos(idxp,:),[1,0,0], 'color','g',  'unit', 'm','alpha',1)
ft_plot_dipole(sourcemodel.pos(idxp,:),[0,1,0], 'color','b',  'unit', 'm','alpha',1)
ft_plot_dipole(sourcemodel.pos(idxp,:),[0,0,1], 'color','r',  'unit', 'm','alpha',1)
view([140 10])

%% Visualization leadfields - dipole leadfields
%--------------------------------------------------------------------------

% sourcemodel_type = 'vol';
sourcemodel_type = 'surf';

switch sourcemodel_type
    case 'vol'
        lf_meg = lf_meg_vol;
        lf_eeg = lf_bem_vol;
%         lf_eeg = lf_fem_vol;
        idxp = idxvol;
    case 'surf'
        lf_meg = lf_meg_surf;
        lf_eeg = lf_bem_surf;
%         lf_eeg = lf_fem_surf;
        idxp = idxsurf;
end


layout_eeg = importdata(fullfile(settings.path2bids,'derivatives',subject,'layout',[subject,'_layout-eeg_orig.mat']));
sensortype = {'mag','planar2','planar3','meg','eeg'};
label      = {'MEG_{MAG}','MEG_{GRAD}','MEG_{GRAD}','MEG','EEG'};
S          = length(sensortype);

[channel] = ft_channelselection('megmag',lf_meg.label);
idx       = contains(lf_meg.label,ft_channelselection('megmag',lf_meg.label));

limits_all    = max(abs(lf_meg.leadfield{idxp}()),[],'all');
idx           = contains(lf_meg.label,ft_channelselection('megmag',lf_meg.label)); % only mags
limits_mag    = max(abs(lf_meg.leadfield{idxp}(idx,:)),[],'all');
idx           = contains(lf_meg.label,ft_channelselection('megplanar',lf_meg.label)); % only grads
limits_planar = max(abs(lf_meg.leadfield{idxp}(idx,:)),[],'all');
limits_eeg    = max(abs(lf_eeg.leadfield{idxp}),[],'all');

% figure('Position', get(0, 'Screensize'))
figure
for didx = 1:3 % loop over direction
        
    for sidx = 1:S % loop over sensortypes
        
        % Generate dummy structure for plotting
        clear dummy
        dummy       = [];
        dummy.time  = 0;
     
        switch sensortype{sidx}
            case 'mag'
                layout      = 'neuromag306mag_helmet.mat';
                coi         = ft_channelselection('megmag',lf_meg.label);
                limits      = limits_mag;
                dummy.avg   = lf_meg.leadfield{idxp}(:,didx);
                dummy.label = lf_meg.label;
                dummy.grad  = grad;
            case 'planar2'
                layout      = 'neuromag306planar_helmet.mat';
                coi         = ft_channelselection('M*2',lf_meg.label);
                limits      = limits_planar;
                dummy.avg   = lf_meg.leadfield{idxp}(:,didx);
                dummy.label = lf_meg.label;
                dummy.grad  = grad;
            case 'planar3'
                layout      = 'neuromag306planar_helmet.mat';
                coi         = ft_channelselection('M*3',lf_meg.label);
                limits      = limits_planar;
                dummy.avg   = lf_meg.leadfield{idxp}(:,didx);
                dummy.label = lf_meg.label;
                dummy.grad  = grad;
            case 'meg'
                layout      = 'neuromag306all_helmet.mat';
                coi         = 'meg';
                limits      = limits_all;
                dummy.avg   = lf_meg.leadfield{idxp}(:,didx);
                dummy.label = lf_meg.label;
                dummy.grad  = grad;
            case 'eeg'
                layout      = layout_eeg;
                coi         = 'all';
                limits      = limits_eeg;
                dummy.avg   = lf_eeg.leadfield{idxp}(:,didx);
                dummy.label = lf_eeg.label;
                dummy.elec  = elec;
        end
        
        cfg             = [];
        cfg.colorbar    = 'yes';
        cfg.zlim        = [-limits,limits]; 
        cfg.comment     = 'no';
        cfg.layout      = layout;
        cfg.interactive = 'no';
        cfg.channel     = coi;
        cfg.figure      = 'gcf'; % plot in current figure
        subplot(3,S,(didx-1)*S+sidx)
        ft_topoplotER(cfg,dummy)
        colormap(flipud(brewermap(64,'RdBu')))
        if didx == 1
            title(label{sidx},'FontSize',15)
        end
    end  
end

%% Visualization leadfields - dipole leadfields (3D)
%--------------------------------------------------------------------------

sourcemodel_type = 'vol';
% sourcemodel_type = 'surf';

switch sourcemodel_type
    case 'vol'
        lf_meg = lf_meg_vol;
        lf_eeg = lf_bem_vol;
%         lf_eeg = lf_fem_vol;
        idxp = idxvol;
    case 'surf'
        lf_meg = lf_meg_surf;
        lf_eeg = lf_bem_surf;
%         lf_eeg = lf_fem_surf;
        idxp = idxsurf;
end

layout_eeg = importdata(fullfile(settings.path2bids,'derivatives',subject,'layout',[subject,'_layout-eeg_orig.mat']));
sensortype = {'mag','planar2','planar3','meg','eeg'};
label      = {'MEG_{MAG}','MEG_{GRAD}','MEG_{GRAD}','MEG','EEG'};
S          = length(sensortype);

% figure('Position', get(0, 'Screensize'))
figure
for didx = 1:3 % loop over direction
        
    for sidx = 1:S % loop over sensortypes
     
        switch sensortype{sidx}
            case 'mag'
                cfg         = [];
                cfg.channel = 'megmag';
                lf          = ft_selectdata(cfg,lf_meg);
                idx         = contains(grad.label,ft_channelselection('megmag',grad.label)); 
            case 'planar2'
                cfg         = [];
                cfg.channel = 'M*2';
                lf          = ft_selectdata(cfg,lf_meg);  
                idx         = contains(grad.label,ft_channelselection('M*2',grad.label)); 
            case 'planar3'
                cfg         = [];
                cfg.channel = 'M*3';
                lf          = ft_selectdata(cfg,lf_meg);   
                idx         = contains(grad.label,ft_channelselection('M*2',grad.label)); 
            case 'meg'
                cfg         = [];
                cfg.channel = 'meg';
                lf          = ft_selectdata(cfg,lf_meg); 
                idx         = contains(grad.label,ft_channelselection('meg',grad.label)); 
            case 'eeg'
                lf = lf_eeg;
        end

        subplot(3,S,(didx-1)*S+sidx)
        switch sensortype{sidx}
            case 'eeg'
                ft_plot_topo3d(elec.chanpos,lf.leadfield{idxp}(:,didx));
            otherwise
                ft_plot_topo3d(grad.chanpos(idx,:),lf.leadfield{idxp}(:,didx));
        end
        colormap(flipud(brewermap(64,'RdBu')))
        if didx == 1
            title(label{sidx},'FontSize',15)
        end

    end
end
       

%% Visualization leadfields - dipole leadfields bem vs. fem
%--------------------------------------------------------------------------

% sourcemodel_type = 'vol';
sourcemodel_type = 'surf';

switch sourcemodel_type
    case 'vol'
        lf_bem = lf_bem_vol;
        lf_fem = lf_fem_vol;
        idxp = idxvol;
    case 'surf'
        lf_bem = lf_bem_surf;
        lf_fem = lf_fem_surf;
        idxp = idxsurf;
end

layout_eeg = importdata(fullfile(settings.path2bids,'derivatives',subject,'layout',[subject,'_layout-eeg_orig.mat']));
headmodels = {'BEM','FEM'};
label      = {'EEG_{BEM}','EEG_{FEM}'};

directions = 3; % all directions
S          = length(headmodels);

limits_bem    =  max(abs(lf_bem.leadfield{idxp}),[],'all');
limits_fem    =  max(abs(lf_fem.leadfield{idxp}),[],'all');

% figure('Position', get(0, 'Screensize'))
figure

% Generate dummy structure for plotting
clear dummy
dummy      = [];
dummy.elec = elec;
dummy.time = 0;
    
for didx = 1:3 % loop over direction
        
    for sidx = 1:2 % eeg headmodels
        
        switch headmodels{sidx}
            case 'BEM'
                limits      = limits_bem;
                dummy.avg   = lf_bem.leadfield{idxp}(:,didx);
                dummy.label = lf_bem.label;
            case 'FEM'
                limits      = limits_fem;
                dummy.avg   = lf_fem.leadfield{idxp}(:,didx);
                dummy.label = lf_fem.label;
        end
            
        cfg             = [];
        cfg.colorbar    = 'yes';
        cfg.zlim        = [-limits,limits]; % commen it out if you don't want limits
        cfg.comment     = 'no';
        cfg.layout      = layout_eeg;
        cfg.interactive = 'no';
        cfg.channel     = 'all';
        cfg.figure      = 'gcf'; % plot in current figure
        subplot(directions,S,(didx-1)*S+sidx)
        ft_topoplotER(cfg,dummy)
        colormap(flipud(brewermap(64,'RdBu')))
        if didx == 1
            title(label{sidx},'FontSize',15)
        end
    end  
end

%% Visualization leadfields - dipole leadfields bem (openmeeg vs. bemcp)
%--------------------------------------------------------------------------

% sourcemodel_type = 'vol';
sourcemodel_type = 'surf';

switch sourcemodel_type
    case 'vol'
        lf_bem   = lf_bem_vol;
        lf_bemcp = lf_bemcp_vol;
        idxp = idxvol;
    case 'surf'
        lf_bem   = lf_bem_surf;
        lf_bemcp = lf_bemcp_surf;
        idxp = idxsurf;
end

layout_eeg = importdata(fullfile(rootpath,'derivatives',subject,'layout',[subject,'_layout-eeg_orig.mat']));
headmodels = {'BEM','BEMCP'};
label      = {'EEG_{BEM}','EEG_{BEMCP}'};

directions = 3; % all directions
S          = length(headmodels);

limits_bem   =  max(abs(lf_bem.leadfield{idxp}),[],'all');
limits_bemcp =  max(abs(lf_bemcp.leadfield{idxp}),[],'all');

figure('Position', get(0, 'Screensize'))

% Generate dummy structure for plotting
clear dummy
dummy      = [];
dummy.elec = elec;
dummy.time = 0;
    
for didx = 1:3 % loop over direction
        
    for sidx = 1:2 % eeg headmodels
        
        switch headmodels{sidx}
            case 'BEM'
                limits      = limits_bem;
                dummy.avg   = lf_bem.leadfield{idxp}(:,didx);
                dummy.label = lf_bem.label;
            case 'BEMCP'
                limits      = limits_bemcp;
                dummy.avg   = lf_bemcp.leadfield{idxp}(:,didx);
                dummy.label = lf_bemcp.label;
        end
            
        cfg             = [];
        cfg.colorbar    = 'yes';
        cfg.zlim        = [-limits,limits]; % commen it out if you don't want limits
        cfg.comment     = 'no';
        cfg.layout      = layout_eeg;
        cfg.interactive = 'no';
        cfg.channel     = 'all';
        cfg.figure      = 'gcf'; % plot in current figure
        subplot(directions,S,(didx-1)*S+sidx)
        ft_topoplotER(cfg,dummy)
        colormap(flipud(brewermap(64,'RdBu')))
        if d == 1
            title(label{sidx},'FontSize',15)
        end
    end  
end




