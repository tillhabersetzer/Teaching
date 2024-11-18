%--------------------------------------------------------------------------
% Generate surface-based sourcemodels of a cortical sheet for mne
%--------------------------------------------------------------------------
% check this tutorial out for reference:
% -> tailored to neuromag coordinate system
% https://www.fieldtriptoolbox.org/workshop/paris2019/handson_anatomy/

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
% define subjects
subjects = [2,3,4];

% resolution - specify only one hemisphere
resolution{1} = '.L.midthickness.4k_fs_LR.surf.gii'; label{1} = '4k';
resolution{2} = '.L.midthickness.8k_fs_LR.surf.gii'; label{2} = '8k';

% option to plot sourcemodel
check = 1;
%--------------------------------------------------------------------------

%% Source model 
%--------------------------------------------------------------------------
for subidx = subjects

    subject   = ['sub-',num2str(subidx,'%02d')];
 
    % Load coordinate transformation matrices
    mri_orig  = ft_read_mri(fullfile(settings.path2bids,'derivatives',subject,'freesurfer',[subject,'.nii']));
    mri_coreg = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_T1w-coregistered.mat'])); % in mm

    trafo_vox2unknown  = mri_orig.transform;
    trafo_vox2neuromag = mri_coreg.transform;
 
    for r = 1:length(resolution)
        
        filename = [subject,resolution{r}];
        path     = fullfile(settings.path2bids,'derivatives',subject,'freesurfer','workbench',filename);
        
        % add right hemisphere
        files       = {path, strrep(path,'.L.','.R.')};
        sourcemodel = ft_read_headshape(files);

        % Coregister mri and sensors
        trafo_unkown2neuromag = trafo_vox2neuromag/trafo_vox2unknown;
        sourcemodel           = ft_transform_geometry(trafo_unkown2neuromag, sourcemodel);
        
        sourcemodel          = ft_determine_units(sourcemodel);
        sourcemodel.coordsys = 'neuromag';
        % sourcemodel.inside   = sourcemodel.atlasroi>0; -> seems not to work with ft_sourcestatistics !
        % sourcemodel          = rmfield(sourcemodel,'atlasroi');
    
        % add inflated surface model for visualization
        inflated          = ft_read_headshape(strrep(files, 'midthickness', 'inflated'));
        inflated          = ft_transform_geometry(trafo_unkown2neuromag, inflated);
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

        % make folder for data
        %---------------------
        dir2save = fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','sourcemodel');
        if ~exist(dir2save, 'dir')
           mkdir(dir2save)
        end

        save(fullfile(dir2save,[subject,'_sourcemodel-corticalsheet',label{r},'.mat']),'sourcemodel'); % in mm
        save(fullfile(dir2save,[subject,'_sourcemodel-inflatedcorticalsheet',label{r},'.mat']),'inflated'); % in mm
     
    end % resolutions
end % subjects      
         