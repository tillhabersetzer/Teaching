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
bp_freq = [4, 120];

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
parfor f = 1:N_files
    path_dataset = [subjectdata.rawdatadir filesep filenames{f} '.fif'];

    %% define trials
    cfg                     = [];
    cfg.dataset             = path_dataset;
    cfg.trialfun            = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = 2;
    cfg.trialdef.prestim    = 0.05;                  % in seconds
    cfg.trialdef.poststim   = 0.29;                   % in seconds
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
    
    %% preprocess data 
    cfg                = [];
    cfg.dataset        = path_dataset;
    cfg.trl            = trl_new;
    cfg.channel        = [{'meg'},strcat('-',subjectdata.badchannels)]; 
    cfg.coilaccuracy   = 0;
    cfg.baselinewindow = [-0.05 0];
    cfg.demean         = 'yes';          
    data               = ft_preprocessing(cfg);
    
    %% reject earlier specified independet components
    if ica_on
        data = reject_independent_components(data,subjectdata,files2preproc);
    end
       
    %% again trial rejection for the rest
    % better to do it after removing of the hugh ecg artifacts
    % magnetometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = {'megmag'}; 
    cfg.artfctdef.threshold.range    = 2500*10^-15; % 1000 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold1]         = ft_artifact_threshold(cfg,data);
    
    % gradiometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = {'megplanar'}; 
    cfg.artfctdef.threshold.range    = 2500*10^-15/(4*10^-2); % 800 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold2]         = ft_artifact_threshold(cfg,data);
    
    cfg                              = [];
    cfg.artfctdef.threshold.artifact = [artifact_threshold1;artifact_threshold2];
    data = ft_rejectartifact(cfg,data);
    
    %% filter data
    cfg           = [];
    cfg.padtype   = 'mirror';
    %cfg.padding    = 1;
    cfg.channel   = 'meg';        
    cfg.bpfilter  = 'yes';
    cfg.bpfreq    = bp_freq;      
    cfg.dftfilter = 'yes';        
    cfg.dftfreq   = [50 100];  
    data          = ft_preprocessing(cfg, data);
    
    %% downsample data
    if downsample_data
        cfg            = [];
        cfg.resamplefs = 400;
        cfg.detrend    = 'no';
        data           = ft_resampledata(cfg,data);
    end
  
    data_preprocessed{f} = data;

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

%% averaging and noise-covariance estimation

channel_right = channel_selection({'MRF','MRT','MRP','MRO'},data_test.label);
channel_left  = channel_selection({'MLF','MLT','MLP','MLO'},data_test.label);
% choose channel
chani = channel_right;
% choose dataset to test
k     = 1;


avg                  = cell(1,N_files);
cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all'; %it will calculate the covariance matrix
                                  % on the timepoints that are before
cfg2         = [];
cfg2.channel = chani;
cfg3         = [];
cfg3.toilim  = [0.06 0.11];
for n = 1:N_files
    data      = ft_selectdata(cfg2,data_preprocessed{n});
    datapst   = ft_redefinetrial(cfg3, data);
    avg{n}    = ft_timelockanalysis(cfg, data);
    avgpst{n} = ft_timelockanalysis(cfg, datapst);
    clear data datapst
end

%% covariance matrix of noise
data_test = avg{k};

cfg                      = [];
cfg.dataset              = [subjectdata.rawdatadir filesep 'empty_pre_tsss.fif'];
cfg.trialfun             = 'ft_trialfun_general'; 
cfg.trialdef.triallength = 10;
cfg.trialdef.ntrials     = inf;
cfg                      = ft_definetrial(cfg);
cfg.channel              = chani;

data_noise               = ft_preprocessing(cfg);

cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg, data_noise);

figure('name','noise covariance')
subplot(1,3,1)
imagesc(avg_noise.cov)
title('empty_pre covariance')
subplot(1,3,2)
imagesc(avg{k}.cov)
title('complete trial covariance')
subplot(1,3,3)
imagesc(avgpst{k}.cov)
title('pst trial covariance')

%% effects of maxfilter
[u,s,v] = svd(avg{k}.cov);
figure;
semilogy(diag(s),'o-');
% get kappa value for beamformer

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

figure
hold on   
ft_plot_headmodel(headmodel_meg,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; %camlight;
alpha 0.4  % make the surface transparent
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); % plot only locations inside the volume
ft_plot_sens(ft_convert_units(data_preprocessed{1}.grad,'mm'),'coilshape', 'circle'); % plot the sensor array
view ([0 -90 0])

%% Compute the leadfield

cfg                 = [];
cfg.channel         = chani;
cfg.grad            = data_test.grad;
cfg.headmodel       = headmodel_meg;
cfg.lcmv.reducerank = 2;                            % default for MEG is 2, for EEG is 3
cfg.sourcemodel     = sourcemodel;
%cfg.normalize      = 'yes';                        % if you are not contrasting
grid                = ft_prepare_leadfield(cfg);

%% first call to sourceanalysis - compute filter
cfg                   = [];
cfg.method            = 'lcmv';
cfg.channel           = chani;
cfg.sourcemodel       = grid;
cfg.headmodel         = headmodel_meg;
cfg.lcmv.keepfilter   = 'yes'; % return filter matrix
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.kappa        = 67;
cfg.lcmv.fixedori     = 'yes';
sourceavg             = ft_sourceanalysis(cfg, avg{k});

%% Perform sourcanalysis
cfg                    = [];
cfg.method             = 'lcmv';
cfg.channel            = chani;
cfg.sourcemodel        = sourcemodel;
cfg.sourcemodel.filter = sourceavg.avg.filter;
cfg.headmodel          = headmodel_meg;
cfg.lcmv.projectnoise  = 'yes';
sourcepst              = ft_sourceanalysis(cfg, avgpst{k});

sourceNAI         = sourcepst;
sourceNAI.avg.pow = sourcepst.avg.pow ./ sourcepst.avg.noise;

%% Interpolate on template mri
templatefile          = 'Z:\MEG_EEG_Toolboxen\fieldtrip-20191127\template\anatomy\single_subj_T1.nii';
template_mri          = ft_read_mri(templatefile);
template_mri.coordsys = 'mni';  % so that FieldTrip knows how to interpret the coordinate system (or 'mni')

sourceavg.pos = template_grid.pos; % alignment with mni space
sourcepst.pos = template_grid.pos; % alignment with mni space
sourceNAI.pos = template_grid.pos; % alignment with mni space
cfg              = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
sourceavg_int    = ft_sourceinterpolate(cfg, sourceavg, template_mri);
sourcepst_int    = ft_sourceinterpolate(cfg, sourcepst, template_mri);
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

[channel] = ft_channelselection('ML',data_preprocessed{1}.label);
channel_right = channel_selection({'MRF-MG','MRT-MG','MRP-MG','MRO-MG'},data_preprocessed{1}.label);
channel_left  = channel_selection({'MLF-MG','MLT-MG','MLP-MG','MLO-MG'},data_preprocessed{1}.label);

%% extra
cfg = [];
for n = 1:N_files
    data_avg{n} = ft_timelockanalysis(cfg,data_preprocessed{n});
end

cfg            = [];
cfg.method     = 'sum';
cfg.updatesens = 'no'; % need old information for dipole fitting
for n = 1:N_files
    data_avg_cmb{n} = ft_combineplanar(cfg,data_avg{n});
end

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
cfg.xlim       = [0, 0.2];
ft_multiplotER(cfg,data_avg_cmb{1});
suptitle('average magnetometer')

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg, data_avg_cmb{1});
suptitle('average combined gradiometers')







