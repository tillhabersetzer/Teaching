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
addpath(fullfile('..','..'))
eval('main_settings')

addpath(fullfile('..','..','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
subidx  = 2;
subject = ['sub-',num2str(subidx,'%02d')];

% Choose timewindow for plotting
% time2plot = [-10 90]; % -10, 90ms
time2plot = [-50 300];

%% Load data 
%--------------------------------------------------------------------------
data = importdata(fullfile(settings.path2project,'derivatives',subject,'sourcelevel',[subject,'_dipolefits-N19mP30m.mat']));

% dipole moments
source_pooled_sym   = data.source_pooled_sym;
source_pooled_nosym = data.source_pooled_nosym;
% vector dipol moments of conditions
source_vec = data.source_vec;
% scalar dipol moments of conditions
source_sca_mean = data.source_sca_mean;
source_sca_svd  = data.source_sca_svd;

conditions      = data.conditions;
chantype        = data.chantype;

sourcemodel   = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_sourcemodel-volumetric.mat']));
sourcemodel   = ft_convert_units(sourcemodel,'m'); 
template_grid = importdata(fullfile(settings.path2fieldtrip,'template','sourcemodel','standard_sourcemodel3d10mm.mat'));
template_grid = ft_convert_units(template_grid,'mm'); 

% positions
%----------
N_chan           = length(chantype);
dippos_sym       = cell(1,N_chan);
dippos_nosym     = cell(1,N_chan);
dippos_sym_mni   = cell(1,N_chan);
dippos_nosym_mni = cell(1,N_chan);

for p=1:N_chan
    dippos_sym{p}   = source_pooled_sym{p}.dip.pos;
    dippos_nosym{p} = source_pooled_nosym{p}.dip.pos;

    % positions in mni space
    idx               = dsearchn(sourcemodel.pos,dippos_sym{p});
    dippos_sym_mni{p} = template_grid.pos(idx,:);

    % positions in mni space
    idx                 = dsearchn(sourcemodel.pos,dippos_nosym{p});
    dippos_nosym_mni{p} = template_grid.pos(idx,:);
end
   
clear data

%% Rescale units - optional for plotting (Am -> nAm)
%--------------------------------------------------------------------------
for cidx = 1:3
    for sidx = 1:N_chan
        source_vec{sidx,cidx}.dip.mom = 10^9*source_vec{sidx,cidx}.dip.mom;
        source_sca_mean{sidx,cidx}    = 10^9*source_sca_mean{sidx,cidx};
        source_sca_svd{sidx,cidx}     = 10^9*source_sca_svd{sidx,cidx};
    end
end

%% (1) Inspect dipole locations
%--------------------------------------------------------------------------

% load mri
%---------
mri = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_T1w-segmented.mat']));
% Keep it in SI-units
mri = ft_convert_units(mri, 'm'); 

for sidx = 1:N_chan
    figure
    hold on
    
    % symmetric dipole fit 
    ft_plot_dipole(source_pooled_sym{sidx}.dip.pos(1,:), mean(source_pooled_sym{sidx}.dip.mom(1:3,:),2), 'color', 'r', 'unit','m')
    ft_plot_dipole(source_pooled_sym{sidx}.dip.pos(2,:), mean(source_pooled_sym{sidx}.dip.mom(4:6,:),2), 'color', 'r', 'unit','m')
    
    % refinement of symmetric dipole fit 
    ft_plot_dipole(source_pooled_nosym{sidx}.dip.pos(1,:), mean(source_pooled_nosym{sidx}.dip.mom(1:3,:),2), 'color', 'g', 'unit','m')
    ft_plot_dipole(source_pooled_nosym{sidx}.dip.pos(2,:), mean(source_pooled_nosym{sidx}.dip.mom(4:6,:),2), 'color', 'g', 'unit','m')
    
    pos = mean(source_pooled_sym{sidx}.dip.pos,1);
    ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.001)
    ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.001)
    ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.001)
    
    ft_plot_crosshair(pos, 'color', [1 1 1]/2);
    
    axis tight
    axis off
    
    view(12, -10)

    subtitle(chantype(sidx))
end

%% (2) Inspect dipole locations in mni space
%--------------------------------------------------------------------------

% this path needs to point to your fieldtrip path
atlas = ft_read_atlas(fullfile(settings.path2fieldtrip,'template','atlas','aal','ROI_MNI_V4.nii')); % mm

cfg            = [];
cfg.atlas      = atlas;
cfg.inputcoord = 'mni';
cfg.output     = 'label';

dippos_loc_sym   = cell(N_chan,2);
dippos_loc_nosym = cell(N_chan,2);

for sidx = 1:N_chan
    for p = 1:2
        positions              = dippos_sym_mni{sidx};
        cfg.roi                = positions(p,:); 
        labels                 = ft_volumelookup(cfg, atlas);
        [~, indx]              = max(labels.count);
        dippos_loc_sym{sidx,p} = labels.name(indx);

        positions                = dippos_nosym_mni{sidx};
        cfg.roi                  = positions(p,:); 
        labels                   = ft_volumelookup(cfg, atlas);
        [~, indx]                = max(labels.count);
        dippos_loc_nosym{sidx,p} = labels.name(indx);
    end
end 

%% (3) Inspect timecourses 
%--------------------------------------------------------------------------

%% All moments, xyz-directions and both hemispheres
%--------------------------------------------------------------------------

for sidx = 1:N_chan
    % Mapping between dipole locations and hemispheres
    %-------------------------------------------------
    pos     = dippos_nosym{sidx};
    % [1,2] if first dipole belongs to left hemisphere, [2,1] if first dipole belongs to right hemisphere
    % 1: left, 2: right
    mapping = check_diploc(pos); 
    
    idxs      = [];
    idxs(1,:) = 1:3; % left
    idxs(2,:) = 4:6; % right
    
    figidxs = [1,3,5;2,4,6];
    name    = {'(left)','(right)'};
    % latency correction
    time    = (source_vec{sidx,1}.time-0.005)*1000; % 5ms
    
    figure
    for i=1:2
        idx    = idxs(mapping(i),:);
        figidx = figidxs(i,:);
        
        mini = min([source_vec{sidx,1}.dip.mom(idx,:),source_vec{sidx,2}.dip.mom(idx,:),source_vec{sidx,3}.dip.mom(idx,:)],[],'all');
        maxi = max([source_vec{sidx,1}.dip.mom(idx,:),source_vec{sidx,2}.dip.mom(idx,:),source_vec{sidx,3}.dip.mom(idx,:)],[],'all');
      
        axisvec = horzcat(time2plot,[mini,maxi]);
        
        subplot(3,2,figidx(1)); 
        plot(time, source_vec{sidx,1}.dip.mom(idx,:), '-')
        if i==1; ylabel('dipole moment / nAm'); end
        legend({'x', 'y', 'z'});
        axis(axisvec) 
        grid on
        title(['click ',name{i}])
        
        subplot(3,2,figidx(2)); 
        plot(time, source_vec{sidx,2}.dip.mom(idx,:), '-')
        if i==1; ylabel('dipole moment / nAm'); end
        legend({'x', 'y', 'z'});
        axis(axisvec)
        grid on
        title(['upchirp ',name{i}])
        
        subplot(3,2,figidx(3)); 
        plot(time, source_vec{sidx,3}.dip.mom(idx,:), '-')
        xlabel('t / ms')
        if i==1; ylabel('dipole moment / nAm'); end
        legend({'x', 'y', 'z'});
        axis(axisvec)
        grid on
        title(['downchirp ',name{i}])
    end
    sgtitle([subject,' | ',chantype{sidx}])
end

%% Visualize the z-component for both hemispheres
%--------------------------------------------------------------------------

for sidx = 1:N_chan
    % Mapping between dipole locations and hemispheres
    %-------------------------------------------------
    pos     = dippos_nosym{sidx};
    mapping = check_diploc(pos); 
    
    idxs = [3,6]; % z-component
    
    name    = {'left hemisphere','right hemisphere'};
    % latency correction
    time    = (source_vec{sidx,1}.time-0.005)*1000; % 5ms
    
    figure
    for i=1:2
    
        idx  = idxs(mapping(i));
        
        mini = min([source_vec{sidx,1}.dip.mom(idx,:),source_vec{sidx,2}.dip.mom(idx,:),source_vec{sidx,3}.dip.mom(idx,:)],[],'all');
        maxi = max([source_vec{sidx,1}.dip.mom(idx,:),source_vec{sidx,2}.dip.mom(idx,:),source_vec{sidx,3}.dip.mom(idx,:)],[],'all');
    
        axisvec = horzcat(time2plot,[mini,maxi]);
        
        subplot(2,1,i); 
        hold on
        plot(time, source_vec{sidx,1}.dip.mom(idx,:), '-')
        plot(time, source_vec{sidx,2}.dip.mom(idx,:), '-')
        plot(time, source_vec{sidx,3}.dip.mom(idx,:), '-')
        legend({'click', 'upchirp', 'downchirp'});
        axis(axisvec)
        grid on
        title(name{i})
        if i==2
            xlabel('t / ms')
        end
        ylabel('dipole moment / nAm')
    end
    sgtitle([subject,': z-component | ',chantype{sidx}])
end

%% Visualize fixed dipole (fixed orientation)
%--------------------------------------------------------------------------
% mean dipolmoment orientation has been used as orientation constraint

for sidx = 1:N_chan
    % Mapping between dipole locations and hemispheres
    %-------------------------------------------------
    pos     = dippos_nosym{sidx};
    mapping = check_diploc(pos); 
    
    name = {'left hemisphere','right hemisphere'};
    % latency correction
    time    = (source_vec{sidx,1}.time-0.005)*1000; % 5ms
    time2plot = [-50,300]; % ms
    
    figure
    for i=1:2
    
        idx  = mapping(i);
        
        mini = min([source_sca_mean{sidx,1}(idx,:),source_sca_mean{sidx,2}(idx,:),source_sca_mean{sidx,3}(idx,:)],[],'all');
        maxi = max([source_sca_mean{sidx,1}(idx,:),source_sca_mean{sidx,2}(idx,:),source_sca_mean{sidx,3}(idx,:)],[],'all');
    
        axisvec = horzcat(time2plot,[mini,maxi]);
        
        subplot(2,1,i); 
        hold on
        plot(time, source_sca_mean{sidx,1}(idx,:), '-')
        plot(time, source_sca_mean{sidx,2}(idx,:), '-')
        plot(time, source_sca_mean{sidx,3}(idx,:), '-')
        legend({'click', 'upchirp', 'downchirp'});
        axis(axisvec)
        grid on
        title(name{i})
        if i ==2
            xlabel('t / ms')
        end
        ylabel('dipole moment / nAm')
    end
    sgtitle([subject,': fixed oriented dipole via mean orientation | ',chantype{sidx}])
end

%% Experimental plots - dont take too serious...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Single timecourse: Projection of moments on direction which explains
%  maximum variance via SVD
%--------------------------------------------------------------------------

for sidx = 1:N_chan
    % Mapping between dipole locations and hemispheres
    %-------------------------------------------------
    pos     = dippos_nosym{sidx};
    mapping = check_diploc(pos); 
    
    idxs      = [];
    idxs(1,:) = 1:3; % left
    idxs(2,:) = 4:6; % right
    
    name = {'left hemisphere','right hemisphere'};
    % latency correction
    time    = (source_vec{sidx,1}.time-0.005)*1000; % 5ms
    
    figure
    for i=1:2
    
        idx = idxs(mapping(i),:);
        
        moments        = source_vec{sidx,1}.dip.mom(idx,:);
        [u, ~, ~]      = svd(moments, 'econ'); 
        moments_clicks = u(:,1)'*moments; 
        
        moments          = source_vec{sidx,2}.dip.mom(idx,:);
        [u, ~, ~]        = svd(moments, 'econ'); 
        moments_upchirps = u(:,1)'*moments; 
        
        moments            = source_vec{sidx,3}.dip.mom(idx,:);
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
    sgtitle([subject,': projection onto maximum variance direction | ',chantype{sidx}])
end

%% Visualize fixed dipole (fixed orientation)
%--------------------------------------------------------------------------
% mean dipolmoment orientation has been used as orientation constraint.
% These orientations were computed via SVD of dipole moment timecourses.
% The orientations explaining largest amount of variance were selected.

%% Visualize fixed dipole (fixed orientation)
%--------------------------------------------------------------------------
% maximum variance dipolmoment orientation via SVD has been used as 
% orientation constraint

for sidx = 1:N_chan
    pos     = dippos_nosym{sidx};
    mapping = check_diploc(pos); 
    
    name = {'left hemisphere','right hemisphere'};
    % latency correction
    time      = (source_vec{sidx,1}.time-0.005)*1000; % 5ms
    
    figure
    for i=1:2
    
        idx  = mapping(i);
        
        mini = min([source_sca_svd{sidx,1}(idx,:),source_sca_svd{sidx,2}(idx,:),source_sca_svd{sidx,3}(idx,:)],[],'all');
        maxi = max([source_sca_svd{sidx,1}(idx,:),source_sca_svd{sidx,2}(idx,:),source_sca_svd{sidx,3}(idx,:)],[],'all');
    
        axisvec = horzcat(time2plot,[mini,maxi]);
        
        subplot(2,1,i); 
        hold on
        plot(time, source_sca_svd{sidx,1}(idx,:), '-')
        plot(time, source_sca_svd{sidx,2}(idx,:), '-')
        plot(time, source_sca_svd{sidx,3}(idx,:), '-')
        legend({'click', 'upchirp', 'downchirp'});
        axis(axisvec)
        grid on
        title(name{i})
        if i ==2
            xlabel('t / ms')
        end
        ylabel('dipole moment / nAm')
    end
    sgtitle([subject,': fixed oriented dipole via svd orientation | ',chantype{sidx}])
end

%% Clean-Up
rmpath(fullfile('..','..'))
rmpath(fullfile('..','..','helper_functions'))

