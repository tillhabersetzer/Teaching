%--------------------------------------------------------------------------
% Generate volumtric sourcemodels for beamforming
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
subjects = [1,2];

% option to plot sourcemodel
check = 1;
%--------------------------------------------------------------------------

%% Source model 
%--------------------------------------------------------------------------

% load template grid
%-------------------
template_grid = importdata(fullfile(settings.path2fieldtrip,'template','sourcemodel','standard_sourcemodel3d10mm.mat'));
template_grid = ft_convert_units(template_grid,'mm'); 

% plot the atlas based grid
if check
    figure
    ft_plot_mesh(template_grid.pos(template_grid.inside,:));
end

% Coregistration of subject specific grid to the atlas based template grid
%-------------------------------------------------------------------------
for subidx = subjects

    subject  = ['sub-',num2str(subidx,'%02d')];

    % make individual subject's grid
    mri_segmented = importdata(fullfile(settings.path2project,'derivatives',subject, ...
                    'forward_modelling',[subject,'_T1w-segmented.mat'])); % mm
    % mri_segmented = ft_convert_units(mri_segmented,'m'); % use SI units
    
    cfg           = [];
    cfg.warpmni   = 'yes';
    cfg.template  = template_grid;
    cfg.nonlinear = 'yes';
    cfg.mri       = mri_segmented;
    cfg.unit      = 'mm';
    sourcemodel   = ft_prepare_sourcemodel(cfg);
    
    % Plot the final source model together with the individual head model and the sensor array
    if check
        % check sourcemodel
        if maxfilter
            megfile = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-clicks-raw_tsss.fif']);
        else
            megfile = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-clicks.fif']);
        end

        shape     = ft_read_headshape(megfile,'unit','mm');
        grad      = ft_convert_units(ft_read_sens(megfile,'senstype','meg'),'mm'); 
        headmodel = importdata(fullfile(settings.path2project,'derivatives',subject, ...
                    'forward_modelling',[subject,'_headmodel-singleshell.mat'])); % mm
        
        figure
        hold on   
        ft_plot_headmodel(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
        alpha 0.5;  % camlight;
        alpha 0.4;  % make the surface transparent
        ft_plot_headshape(shape);
        ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); % plot only locations inside the volume
        ft_plot_sens(grad, 'style', '*b');
        view ([0 -90 0])
        title(subject)
    end
    
    % Save data
    %----------   
    % make folder for data
    dir2save = fullfile(settings.path2project,'derivatives',subject,'forward_modelling');
    if ~exist(dir2save, 'dir')
       mkdir(dir2save)
    end
    save(fullfile(dir2save,[subject,'_sourcemodel-volumetric.mat']),'sourcemodel'); % in mm

end 

%% Clean-Up
rmpath(fullfile('..','..'))
