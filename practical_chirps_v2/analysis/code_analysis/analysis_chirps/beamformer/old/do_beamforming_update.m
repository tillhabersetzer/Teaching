close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the avarages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject02'};

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [0.5, 150];

% apply ica
ica_on = 1;

% downsample data
downsample_data = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);
% addpath for ica functions
addpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);

%% initialize

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj
    
% subject selected
subject = subjects{s};
eval(subject)

filenames       = get_filenames(subjectdata,files2preproc);
F               = length(filenames);
filename4saving = [filenames,{files2preproc}];
%     path2save = subjectdata.ica_stories;
N_files         = F;

% loop over selected files
%-------------------------
data_preprocessed = cell(1,N_files);

%parfor f = 1:N_files
for f = 1:N_files
    path_dataset = [subjectdata.rawdatadir filesep filenames{f} '.fif'];
    %% preprocess data
    cfg              = [];
    cfg.dataset      = path_dataset;
    % channels of interest and exclude bad channels
    cfg.channel      = [{'meg'},strcat('-',subjectdata.badchannels)]; 
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;
    data_cont        = ft_preprocessing(cfg);
    
    %% reject earlier specified independet components
    if ica_on
        data_cont = reject_independent_components(data_cont,subjectdata,files2preproc);
    end
    
    %% filter data    
    cfg            = [];
    cfg.channel    = 'meg';        
    cfg.detrend    = 'yes';
    cfg.demean     = 'yes';         
    cfg.continuous = 'yes';
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = bp_freq;      
    cfg.dftfilter  = 'yes';        
    cfg.dftfreq    = [50 100 150];  
    data_filt      = ft_preprocessing(cfg, data_cont);
    
    clear data_cont

    %% define trials
    cfg                     = [];
    cfg.dataset             = path_dataset;
    cfg.trialfun            = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = 2;
    cfg.trialdef.prestim    = 0.05;                  % in seconds
    cfg.trialdef.poststim   = 0.3;                   % in seconds
    cfg                     = ft_definetrial(cfg);
    trl                     = cfg.trl;
    
    %% detect bad trials
    cfg                        = [];
    cfg.trl                    = trl;
    cfg.dataset                = path_dataset;
    cfg.artfctdef.jump.channel = [{'meg'},strcat('-',subjectdata.badchannels)]; 
    [~, artifact_jump]         = ft_artifact_jump(cfg);
    
    cfg                        = [];
    cfg.trl                    = trl;
    cfg.dataset                = path_dataset;
    cfg.artfctdef.clip.channel = [{'meg'},strcat('-',subjectdata.badchannels)]; 
    [~, artifact_clip]         = ft_artifact_clip(cfg);
    
    cfg                         = [];
    cfg.trl                     = trl;
    cfg.dataset                 = path_dataset;
    cfg.artfctdef.jump.artifact = artifact_jump;
    cfg.artfctdef.clip.artifact = artifact_clip;
    cfg = ft_rejectartifact(cfg);
    trl_new = cfg.trl;
     
    %% apply trial definition
    cfg     = [];
    cfg.trl = trl_new;
    data    = ft_redefinetrial(cfg, data_filt);
    
    clear data_filt
    
    %% again trial rejection for the rest
    % better to do it after removing of the hugh ecg artifacts
    % magnetometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = {'megmag'}; 
    cfg.artfctdef.threshold.range    = 2000*10^-15; % 1000 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold1]         = ft_artifact_threshold(cfg,data);
    
    % gradiometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = {'megplanar'}; 
    cfg.artfctdef.threshold.range    = 2000*10^-15/(4*10^-2); % 800 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold2]         = ft_artifact_threshold(cfg,data);
    
    cfg                              = [];
    cfg.artfctdef.threshold.artifact = [artifact_threshold1;artifact_threshold2];
    data = ft_rejectartifact(cfg,data);
    
    %% downsample data
    if downsample_data
        cfg            = [];
        cfg.resamplefs = 400;
        cfg.detrend    = 'no';
        data           = ft_resampledata(cfg,data);
    end
  
    data_preprocessed{f} = data;
    
    clear data

end % files

%% combine data
cfg                = [];
cfg.keepsampleinfo = 'no';
data_all           = ft_appenddata(cfg,data_preprocessed{:});

% average sensor array
sens          = cellfun(@(v) v.grad,data_preprocessed);
[asens, afid] = ft_average_sens(sens);
data_all.grad = asens;

% add combined file
data_preprocessed{N_files+1} = data_all;
N_files                      = length(data_preprocessed);

%% only look at gradiometers
data_megplanar  = cell(size(data_preprocessed));
cfg         = [];
cfg.channel = 'megplanar';
for n = 1:N_files
    data_megplanar{n} = ft_selectdata(cfg,data_preprocessed{n});
end
%% new try - sourcemodel
template_grid = importdata(['Z:\MEG_EEG_Toolboxen\fieldtrip-20191127\' ...
                            'template\sourcemodel\standard_sourcemodel3d10mm.mat']);
template_grid = ft_convert_units(template_grid,'mm');

%% load atlas and create binary mask (restrict your template grid further with a mask)
atlas = ft_read_atlas('Z:\MEG_EEG_Toolboxen\fieldtrip-20191127\template\atlas\aal\ROI_MNI_V4.nii');

cfg            = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;
cfg.inputcoord = 'mni';
mask           = ft_volumelookup(cfg,template_grid);

% create temporary mask according to the atlas entries
tmp         = repmat(template_grid.inside,1,1);
tmp(tmp==1) = 0;
tmp(mask)   = 1;

% define inside locations according to the atlas based mask
template_grid.inside = tmp;

% plot the atlas based grid
figure; 
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%% Inverse-warp the subject specific grid to the atlas based template grid (make individual subject's grid)
mri_segmented = importdata([subjectdata.headmodel filesep 'mri_segmented.mat']);
%mri_segmented_cm = ft_convert_units(mri_segmented, 'cm');

cfg           = [];
cfg.warpmni   = 'yes';
cfg.template  = template_grid;
cfg.nonlinear = 'yes';
cfg.mri       = mri_segmented;
cfg.unit      = 'mm';
sourcemodel   = ft_prepare_sourcemodel(cfg);

%% Plot the final source model together with the individual head model and the sensor array
headmodel_meg = importdata([subjectdata.headmodel filesep 'headmodel_meg.mat']);
%headmodel_meg = ft_convert_units(headmodel_meg, 'm');
%sourcemodel   = ft_convert_units(sourcemodel, 'm');

figure
hold on   
ft_plot_headmodel(headmodel_meg,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; %camlight;
alpha 0.4  % make the surface transparent
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); % plot only locations inside the volume
ft_plot_sens(ft_convert_units(data_preprocessed{1}.grad,'mm'),'coilshape', 'circle'); % plot the sensor array
view ([0 -90 0])

%% Compute the leadfield
data_test          = data_pre
cfg                = [];
% cfg.channel        = data_preprocessed{1}.label;   % ensure that rejected sensors are not present
% cfg.grad           = data_preprocessed{1}.grad;
cfg.channel        = data_test.label;   % ensure that rejected sensors are not present
cfg.grad           = data_test.grad;
cfg.headmodel      = headmodel_meg;
cfg.sam.reducerank = 2;                            % default for MEG is 2, for EEG is 3
cfg.sourcemodel    = sourcemodel;
%cfg.normalize      = 'yes';                        % if you are not contrasting
grid               = ft_prepare_leadfield(cfg);

%% compute data covariance
% around M100 component
% 
% cfg        = [];
% cfg.toilim = [-0.05 0.0];
% datapre{1} = ft_redefinetrial(cfg, data_preprocessed{1});

cfg         = [];
cfg.toilim  = [0.05 0.18];
datapost{5} = ft_redefinetrial(cfg, data_megplanar{5});

cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg{5}               = ft_timelockanalysis(cfg,data_megplanar{5});

cfg            = [];
cfg.covariance = 'yes';
% avgpre{1}      = ft_timelockanalysis(cfg,datapre{1});
avgpst{5}      = ft_timelockanalysis(cfg,datapost{5});

%% first call to sourceanalysis - compute filter
cfg               = [];
cfg.method        = 'lcmv';
%cfg.channel      = ft_channelselection('megplanar',data_preprocessed{1}.label);
cfg.sourcemodel   = grid;
cfg.headmodel     = headmodel_meg;
cfg.keepfilter    = 'yes'; % return filter matrix
cfg.projectnoise  = 'yes';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.kappa    = 60;  % rank of data
sourceavg{5}      = ft_sourceanalysis(cfg, avg{5});

%% Perform sourcanalysis
cfg                    = [];
cfg.method             = 'lcmv';
%cfg.channel            = ft_channelselection('megplanar',data_preprocessed{1}.label);
cfg.sourcemodel        = sourcemodel;
cfg.sourcemodel.filter = sourceavg{5}.avg.filter;
cfg.headmodel          = headmodel_meg;
%sourcepreS1{1}         = ft_sourceanalysis(cfg, avgpre{1});
sourcepst{5}           = ft_sourceanalysis(cfg, avgpst{5});

sourceNAI         = sourceavg{5};
sourceNAI.avg.pow = sourceavg{5}.avg.pow ./ sourceavg{5}.avg.noise;

%% Subtract the two conditions
% cfg           = [];
% cfg.parameter = 'avg.pow';
% cfg.operation = '((x1-x2)./x2)*100';  % data is expressed in percentage change from pre stimulus baseline
% contrast{1}   = ft_math(cfg,sourcepstS1{1},sourcepreS1{1});

%% Interpolate on template mri
templatefile          = 'Z:\MEG_EEG_Toolboxen\fieldtrip-20191127\template\anatomy\single_subj_T1.nii';
template_mri          = ft_read_mri(templatefile);
template_mri.coordsys = 'mni';  % so that FieldTrip knows how to interpret the coordinate system (or 'mni')

sourceavg{5}.pos = template_grid.pos; % alignment with mni space
sourcepst{5}.pos = template_grid.pos; % alignment with mni space
sourceNAI.pos    = template_grid.pos; % alignment with mni space
cfg              = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
sourceavg_int    = ft_sourceinterpolate(cfg, sourceavg{5}, template_mri);
sourcepst_int    = ft_sourceinterpolate(cfg, sourcepst{5}, template_mri);
sourceNAI_int    = ft_sourceinterpolate(cfg, sourceNAI, template_mri);
%% Plot result
cfg              = [];
cfg.atlas        = atlas;
cfg.method       = 'ortho';
cfg.funparameter = 'pow';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.latency      = [0.07,0.1];
ft_sourceplot(cfg,sourceavg_int);
ft_sourceplot(cfg,sourcepst_int);
ft_sourceplot(cfg,sourceNAI_int);

%% Plot the result in parceled brain space

cfg           = [];
cfg.parameter = 'tissue';
mri2          = ft_sourceinterpolate(cfg, atlas, template_mri);










end