%--------------------------------------------------------------------------
% (1) Inspect dipole fits
% (2) Inspect dipole locations in mni space
% (3) Inspect dipole timecourse
%   - all moments xyz directions
%   - z directions only
%   - fixed oriented dipole (like Stefan did in paper?!)
%   Experimental plots
%   - projection onto maximum variance orientation
%   - fixed oriented dipole with orientation via SVD (recommendation)
%--------------------------------------------------------------------------

close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','..','subjectdata'))
eval('main_settings')

addpath(fullfile('..','..','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
subidx  = 3;
subject = ['sub-',num2str(subidx,'%02d')];

% choose channel type for plotting
chantype = 'meg';
% chantype = 'eeg';

% Choose which type of evoked field (erf) to fit
erf_type = 'N19mP30m'; % type for practical
% erf_type = 'N100m'; 

% Choose timewindow for plotting
% time2plot = [-10 90]; % -10, 90ms
time2plot = [-50 300];

%% Load data 
%--------------------------------------------------------------------------
data = importdata(fullfile(settings.path2bids,'derivatives',subject,'sourcelevel',[subject,'_dipolefits-',erf_type,'_',chantype,'.mat']));

% dipole moments
source_pooled       = data.source_pooled;
source_pooled_nosym = data.source_pooled_nosym;
% vector dipol moments of conditions
source_vec = data.source_vec;
% scalar dipol moments of conditions
source_sca_mean = data.source_sca_mean;
source_sca_svd  = data.source_sca_svd;

conditions       = data.conditions;
sourcemodel_type = data.sourcemodel_type;

% positions
%----------
dippos = data.dippos;
switch sourcemodel_type
    case 'vol'
        sourcemodel   = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','sourcemodel',[subject,'_sourcemodel-volumetric.mat']));
        template_grid = importdata(fullfile(settings.path2fieldtrip,'template','sourcemodel','standard_sourcemodel3d10mm.mat'));
        template_grid = ft_convert_units(template_grid,'mm'); 
        % positions in mni space
        idx         = dsearchn(sourcemodel.pos,data.dippos.pos_sym);
        pos_sym_mni = template_grid.pos(idx,:);
        
        % positions in mni space
        idx           = dsearchn(sourcemodel.pos,data.dippos.pos_nosym);
        pos_nosym_mni = template_grid.pos(idx,:);
    case 'surf'
        sourcemodel = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','sourcemodel',[subject,'_sourcemodel-corticalsheet4k.mat']));
end


clear data

%% Rescale units - optional for plotting (Am -> nAm)
%--------------------------------------------------------------------------
for cidx = 1:3
    source_vec{cidx}.dip.mom = 10^9*source_vec{cidx}.dip.mom;
    source_sca_mean{cidx}    = 10^9*source_sca_mean{cidx};
    source_sca_svd{cidx}     = 10^9*source_sca_svd{cidx};
end

%% (1) Inspect dipole locations
%--------------------------------------------------------------------------

% load mri
%---------
mri = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_T1w-segmented.mat']));
% Keep it in SI-units
mri = ft_convert_units(mri, 'm'); 

figure
hold on

% symmetric dipole fit 
ft_plot_dipole(source_pooled.dip.pos(1,:), mean(source_pooled.dip.mom(1:3,:),2), 'color', 'r', 'unit','m')
ft_plot_dipole(source_pooled.dip.pos(2,:), mean(source_pooled.dip.mom(4:6,:),2), 'color', 'r', 'unit','m')

% refinement of symmetric dipole fit 
ft_plot_dipole(source_pooled_nosym.dip.pos(1,:), mean(source_pooled_nosym.dip.mom(1:3,:),2), 'color', 'g', 'unit','m')
ft_plot_dipole(source_pooled_nosym.dip.pos(2,:), mean(source_pooled_nosym.dip.mom(4:6,:),2), 'color', 'g', 'unit','m')

pos = mean(source_pooled.dip.pos,1);
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.001)
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.001)
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.001)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis tight
axis off

view(12, -10)

%% (2) Inspect dipole locations in mni space
%--------------------------------------------------------------------------

switch sourcemodel_type
    case 'vol'
        % this path needs to point to your fieldtrip path
        atlas = ft_read_atlas(fullfile(settings.path2fieldtrip,'template','atlas','aal','ROI_MNI_V4.nii')); % mm
        
        cfg            = [];
        cfg.atlas      = atlas;
        cfg.inputcoord = 'mni';
        cfg.output     = 'label';
        
        positions = vertcat(pos_sym_mni,pos_nosym_mni);
        dip_loc   = cell(4,1);
        for p=1:4
            cfg.roi    = positions(p,:); 
            labels     = ft_volumelookup(cfg, atlas);
            [~, indx]  = max(labels.count);
            dip_loc{p} = labels.name(indx);
        end
    case 'surf'
        atlasname  = '.L.aparc.a2009s.4k_fs_LR.label.gii';
        atlas      = generate_surfatlas(atlasname,sourcemodel,settings.path2bids,subject);

        positions = vertcat(dippos.pos_sym,dippos.pos_nosym)*1000; % m -> mm
        dip_loc   = cell(4,1);
        for p=1:4
            idx        = dsearchn(atlas.pos,positions(p,:));
            dip_loc{p} = atlas.parcellationlabel(atlas.parcellation(idx));
        end
end

%% (3) Inspect timecourses 
%--------------------------------------------------------------------------

%% All moments, xyz-directions and both hemispheres
%--------------------------------------------------------------------------

% Mapping between dipole locations and hemispheres
%-------------------------------------------------
pos     = dippos.pos_nosym;
% [1,2] if first dipole belongs to left hemisphere, [2,1] if first dipole belongs to right hemisphere
% 1: left, 2: right
mapping = check_diploc(pos); 

idxs      = [];
idxs(1,:) = 1:3; % left
idxs(2,:) = 4:6; % right

figidxs = [1,3,5;2,4,6];
name    = {'(left)','(right)'};
% latency correction
time    = (source_vec{1}.time-0.005)*1000; % 5ms

figure
for i=1:2
    idx    = idxs(mapping(i),:);
    figidx = figidxs(i,:);
    
    mini = min([source_vec{1}.dip.mom(idx,:),source_vec{2}.dip.mom(idx,:),source_vec{3}.dip.mom(idx,:)],[],'all');
    maxi = max([source_vec{1}.dip.mom(idx,:),source_vec{2}.dip.mom(idx,:),source_vec{3}.dip.mom(idx,:)],[],'all');
  
    axisvec = horzcat(time2plot,[mini,maxi]);
    
    subplot(3,2,figidx(1)); 
    plot(time, source_vec{1}.dip.mom(idx,:), '-')
    if i==1; ylabel('dipole moment / nAm'); end
    legend({'x', 'y', 'z'});
    axis(axisvec) 
    grid on
    title(['click ',name{i}])
    
    subplot(3,2,figidx(2)); 
    plot(time, source_vec{2}.dip.mom(idx,:), '-')
    if i==1; ylabel('dipole moment / nAm'); end
    legend({'x', 'y', 'z'});
    axis(axisvec)
    grid on
    title(['upchirp ',name{i}])
    
    subplot(3,2,figidx(3)); 
    plot(time, source_vec{3}.dip.mom(idx,:), '-')
    xlabel('t / ms')
    if i==1; ylabel('dipole moment / nAm'); end
    legend({'x', 'y', 'z'});
    axis(axisvec)
    grid on
    title(['downchirp ',name{i}])
end
sgtitle(subject)

%% Visualize the z-component for both hemispheres
%--------------------------------------------------------------------------

% Mapping between dipole locations and hemispheres
%-------------------------------------------------
pos     = dippos.pos_nosym;
mapping = check_diploc(pos); 

idxs = [3,6]; % z-component

name    = {'left hemisphere','right hemisphere'};
% latency correction
time    = (source_vec{1}.time-0.005)*1000; % 5ms

figure
for i=1:2

    idx  = idxs(mapping(i));
    
    mini = min([source_vec{1}.dip.mom(idx,:),source_vec{2}.dip.mom(idx,:),source_vec{3}.dip.mom(idx,:)],[],'all');
    maxi = max([source_vec{1}.dip.mom(idx,:),source_vec{2}.dip.mom(idx,:),source_vec{3}.dip.mom(idx,:)],[],'all');

    axisvec = horzcat(time2plot,[mini,maxi]);
    
    subplot(2,1,i); 
    hold on
    plot(time, source_vec{1}.dip.mom(idx,:), '-')
    plot(time, source_vec{2}.dip.mom(idx,:), '-')
    plot(time, source_vec{3}.dip.mom(idx,:), '-')
    legend({'click', 'upchirp', 'downchirp'});
    axis(axisvec)
    grid on
    title(name{i})
    if i==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm')
end
sgtitle([subject,': z-component'])

%% Visualize fixed dipole (fixed orientation)
%--------------------------------------------------------------------------
% mean dipolmoment orientation has been used as orientation constraint

% Mapping between dipole locations and hemispheres
%-------------------------------------------------
pos     = dippos.pos_nosym;
mapping = check_diploc(pos); 

name = {'left hemisphere','right hemisphere'};
% latency correction
time    = (source_vec{1}.time-0.005)*1000; % 5ms
time2plot = [-50,300]; % ms

figure
for i=1:2

    idx  = mapping(i);
    
    mini = min([source_sca_mean{1}(idx,:),source_sca_mean{2}(idx,:),source_sca_mean{3}(idx,:)],[],'all');
    maxi = max([source_sca_mean{1}(idx,:),source_sca_mean{2}(idx,:),source_sca_mean{3}(idx,:)],[],'all');

    axisvec = horzcat(time2plot,[mini,maxi]);
    
    subplot(2,1,i); 
    hold on
    plot(time, source_sca_mean{1}(idx,:), '-')
    plot(time, source_sca_mean{2}(idx,:), '-')
    plot(time, source_sca_mean{3}(idx,:), '-')
    legend({'click', 'upchirp', 'downchirp'});
    axis(axisvec)
    grid on
    title(name{i})
    if i ==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm')
end
sgtitle([subject,': fixed oriented dipole via mean orientation'])

%% Experimental plots - dont take too serious...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Single timecourse: Projection of moments on direction which explains
%  maximum variance via SVD
%--------------------------------------------------------------------------

% Mapping between dipole locations and hemispheres
%-------------------------------------------------
pos     = dippos.pos_nosym;
mapping = check_diploc(pos); 

idxs      = [];
idxs(1,:) = 1:3; % left
idxs(2,:) = 4:6; % right

name = {'left hemisphere','right hemisphere'};
% latency correction
time    = (source_vec{1}.time-0.005)*1000; % 5ms

figure
for i=1:2

    idx = idxs(mapping(i),:);
    
    moments        = source_vec{1}.dip.mom(idx,:);
    [u, ~, ~]      = svd(moments, 'econ'); 
    moments_clicks = u(:,1)'*moments; 
    
    moments          = source_vec{2}.dip.mom(idx,:);
    [u, ~, ~]        = svd(moments, 'econ'); 
    moments_upchirps = u(:,1)'*moments; 
    
    moments            = source_vec{3}.dip.mom(idx,:);
    [u, ~, ~]          = svd(moments, 'econ'); 
    moments_downchirps = u(:,1)'*moments; 

    mini = min([moments_clicks,moments_upchirps,moments_downchirps],[],'all');
    maxi = max([moments_clicks,moments_upchirps,moments_downchirps],[],'all');

    axisvec = horzcat(time2plot,[mini,maxi]);
    
    subplot(2,1,i)
    hold on
    plot(time, moments_clicks, '-')
    plot(time, moments_upchirps, '-')
    plot(time, moments_downchirps, '-')
    legend({'click', 'upchirp', 'downchirp'});
    axis(axisvec)
    grid on
    title(name{i})
    if i ==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm')
end
sgtitle([subject,': projection onto maximum variance direction'])

%% Visualize fixed dipole (fixed orientation)
%--------------------------------------------------------------------------
% mean dipolmoment orientation has been used as orientation constraint.
% These orientations were computed via SVD of dipole moment timecourses.
% The orientations explaining largest amount of variance were selected.

%% Visualize fixed dipole (fixed orientation)
%--------------------------------------------------------------------------
% maximum variance dipolmoment orientation via SVD has been used as 
% orientation constraint

pos     = dippos.pos_nosym;
mapping = check_diploc(pos); 

name = {'left hemisphere','right hemisphere'};
% latency correction
time      = (source_vec{1}.time-0.005)*1000; % 5ms


figure
for i=1:2

    idx  = mapping(i);
    
    mini = min([source_sca_svd{1}(idx,:),source_sca_svd{2}(idx,:),source_sca_svd{3}(idx,:)],[],'all');
    maxi = max([source_sca_svd{1}(idx,:),source_sca_svd{2}(idx,:),source_sca_svd{3}(idx,:)],[],'all');

    axisvec = horzcat(time2plot,[mini,maxi]);
    
    subplot(2,1,i); 
    hold on
    plot(time, source_sca_svd{1}(idx,:), '-')
    plot(time, source_sca_svd{2}(idx,:), '-')
    plot(time, source_sca_svd{3}(idx,:), '-')
    legend({'click', 'upchirp', 'downchirp'});
    axis(axisvec)
    grid on
    title(name{i})
    if i ==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm')
end
sgtitle([subject,': fixed oriented dipole via svd orientation'])

