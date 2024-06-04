 close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the averages and
% beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject04'};

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [4, 30];

% apply ica
ica_on = 0;

% choose files for detected ica components
ica_files = 'all_maxfilter';

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

filenames = get_filenames(subjectdata,files2preproc);
N_files   = length(filenames);

% loop over selected files
%-------------------------
data_preprocessed = cell(1,N_files);

%parfor f = 1:N_files
for f = 1:N_files
    
    path_dataset = [subjectdata.rawdatadir filesep filenames{f} '.fif'];
    % filter continuous data to avoid edge artifacts
    cfg              = [];
    cfg.dataset      = path_dataset;
    cfg.channel      = 'meg'; 
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;            % ensure that sensors are expressed in SI units
    data             = ft_preprocessing(cfg);   
    
    %% reject earlier specified independet components
    if ica_on
        data = reject_independent_components(data,subjectdata,ica_files);
    end
    
    %% filter data
    cfg              = [];
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = bp_freq;
    cfg.dftfilter    = 'yes';        % enable notch filtering to eliminate power line noise
    cfg.dftfreq      = [50 100 150]; % set up the frequencies for notch filtering
    cfg.coilaccuracy = 0;
    data             = ft_preprocessing(cfg,data);   

    %% define trials
    cfg                     = [];
    cfg.dataset             = path_dataset;
    cfg.trialfun            = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = 2;
    cfg.trialdef.prestim    = 0.5;                  % in seconds
    cfg.trialdef.poststim   = 0.5;                  % in seconds
    cfg                     = ft_definetrial(cfg);
    trl                     = cfg.trl;
    
    %% detect bad trials
    cfg                        = [];
    cfg.trl                    = trl;
    cfg.dataset                = path_dataset;
    cfg.artfctdef.jump.channel = 'meg'; 
    [~, artifact_jump]         = ft_artifact_jump(cfg);
    
    cfg                        = [];
    cfg.trl                    = trl;
    cfg.dataset                = path_dataset;
    cfg.artfctdef.clip.channel = 'meg'; 
    [~, artifact_clip]         = ft_artifact_clip(cfg);
    
    cfg                         = [];
    cfg.trl                     = trl;
    cfg.dataset                 = path_dataset;
    cfg.artfctdef.jump.artifact = artifact_jump;
    cfg.artfctdef.clip.artifact = artifact_clip;
    cfg = ft_rejectartifact(cfg);
    trl_new = cfg.trl;
    
    %% epoch data
    cfg                = [];
    cfg.trl            = trl_new;          
    data               = ft_redefinetrial(cfg,data);
    
    %% baseline correction
%     cfg                = [];
%     cfg.baselinewindow = [-0.5 0];
%     cfg.demean         = 'yes';
%     data               = ft_preprocessing(cfg,data); 
    
    %% again trial rejection for the rest
    % better to do it after removing of the hugh ecg artifacts
    % magnetometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = 'megmag'; 
    cfg.artfctdef.threshold.range    = 2500*10^-15; % 1000 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold1]         = ft_artifact_threshold(cfg,data);
    
    % gradiometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = 'megplanar'; 
    cfg.artfctdef.threshold.range    = 2500*10^-15/(4*10^-2); % 800 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold2]         = ft_artifact_threshold(cfg,data);
    
    cfg                              = [];
    cfg.artfctdef.threshold.artifact = [artifact_threshold1;artifact_threshold2];
    data = ft_rejectartifact(cfg,data);
    
    %% downsample data
    if downsample_data
        cfg            = [];
        cfg.resamplefs = 300;
        cfg.detrend    = 'no';
        data           = ft_resampledata(cfg,data);
    end
  
    data_preprocessed{f} = data;
    
    clear data

end % files
    
%% append data
hdr                   = data_preprocessed{1}.hdr;
cfg                   = [];
cfg.keepsampleinfo    = 'no';
data_preprocessed     = ft_appenddata(cfg,data_preprocessed{:});
data_preprocessed.hdr = hdr;

%% sourcemodel

% sourcemodel
template_grid = importdata(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ...
                            'template\sourcemodel\standard_sourcemodel3d10mm.mat']);
template_grid = ft_convert_units(template_grid,'mm'); 
% template_grid = ft_convert_units(template_grid,'m'); % use SI units

% load atlas and create binary mask (restrict your template grid further with a mask)
atlas = ft_read_atlas('Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\template\atlas\aal\ROI_MNI_V4.nii');
% atlas = ft_convert_units(atlas,'m'); % use SI units

cfg            = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;
cfg.inputcoord = 'mni';
mask           = ft_volumelookup(cfg,template_grid);

% create temporary mask according to the atlas entries
tmp         = repmat(template_grid.inside,1,1);
tmp(tmp==1) = 0;
tmp(mask)   = 1;

% or use instead - seems not to work!
% template_grid.inside = false(template_grid.dim);
% template_grid.inside(mask==1) = true;

% define inside locations according to the atlas based mask
template_grid.inside = tmp;

% plot the atlas based grid
figure; 
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

% Inverse-warp the subject specific grid to the atlas based template grid (make individual subject's grid)
mri_segmented = importdata([subjectdata.headmodel filesep 'mri_segmented.mat']);% mm
% mri_segmented = ft_convert_units(mri_segmented,'m'); % use SI units

cfg           = [];
cfg.warpmni   = 'yes';
cfg.template  = template_grid;
cfg.nonlinear = 'yes';
cfg.mri       = mri_segmented;
cfg.unit      = 'mm';
sourcemodel   = ft_prepare_sourcemodel(cfg);

% Plot the final source model together with the individual head model and the sensor array
headmodel_meg = importdata([subjectdata.headmodel filesep 'headmodel_meg.mat']); % mm

figure
hold on   
ft_plot_headmodel(headmodel_meg,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; %camlight;
alpha 0.4;  % make the surface transparent
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); % plot only locations inside the volume
ft_plot_sens(ft_convert_units(data_preprocessed.grad,'mm'),'coilshape', 'circle'); % plot the sensor array
view ([0 -90 0])

% Compute the leadfield
% data_preprocessed_mm = ft_convert_units(data_preprocessed,'mm'); % does not work

%% try new stuff
%--------------
cfg                = [];
cfg.channel        = 'megplanar';
data_preprocessed1 = ft_selectdata(cfg,data_preprocessed);

cfg             = [];
cfg.channel     = 'megplanar';
cfg.grad        = data_preprocessed1.grad;
cfg.headmodel   = headmodel_meg;
cfg.reducerank  = 2;                            % default for MEG is 2, for EEG is 3
cfg.sourcemodel = sourcemodel;
%cfg.normalize   = 'yes';                       % if you are not contrasting
grid            = ft_prepare_leadfield(cfg);

cfg        = [];
%cfg.toilim = [-.18 -.05];
cfg.toilim = [-0.5 -0.2];
datapre    = ft_redefinetrial(cfg, data_preprocessed1);
%cfg.toilim = [.05 .18];
cfg.toilim = [0 0.3];
datapost   = ft_redefinetrial(cfg, data_preprocessed1);

cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = [-0.5 0.5];
%cfg.covariancewindow = [0 0.5];
avg                  = ft_timelockanalysis(cfg,data_preprocessed1);

% figure
% imagesc(avg.cov)

cfg            = [];
cfg.covariance = 'yes';
avgpre         = ft_timelockanalysis(cfg,datapre);
avgpst         = ft_timelockanalysis(cfg,datapost);

% spatial filter on the basis of entire data
cfg             = [];
cfg.method      = 'lcmv';
%cfg.method      = 'sam';
cfg.sourcemodel = grid;
cfg.headmodel   = headmodel_meg;
cfg.lcmv.keepfilter  = 'yes';
%cfg.sam.keepfilter  = 'yes';
%cfg.channel    = 'megplanar';
sourceavg       = ft_sourceanalysis(cfg, avg);


% reconstruct activity with precomputed filter
cfg                    = [];
cfg.method             = 'lcmv';
%cfg.method             = 'sam';
cfg.sourcemodel        = grid;
cfg.sourcemodel.filter = sourceavg.avg.filter;
cfg.headmodel          = headmodel_meg;
sourcepreS1            = ft_sourceanalysis(cfg, avgpre);
sourcepstS1            = ft_sourceanalysis(cfg, avgpst);

cfg           = [];
cfg.parameter = 'avg.pow';
cfg.operation = '((x1-x2)./x2)*100'; %percentage change from baseline
S1bl          = ft_math(cfg,sourcepstS1,sourcepreS1);

templatefile     = 'Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\template\anatomy\single_subj_T1.nii';
template_mri     = ft_read_mri(templatefile);
S1bl.pos         = template_grid.pos; % template grid position are in mni space
%sourcepstS1.pos  = template_grid.pos; % template grid position are in mni space
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int       = ft_sourceinterpolate(cfg, S1bl, template_mri);
%source_int       = ft_sourceinterpolate(cfg, sourcepstS1, template_mri);

% look at plot
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.location      = [64 -32 8];
cfg.funcolormap   = 'jet';
ft_sourceplot(cfg,source_int);

%% other new stuff
templatefile     = 'Z:/Software/MEG_EEG_Toolboxen/fieldtrip-20191127/external/spm8/templates/T1.nii';
template_mri     = ft_read_mri(templatefile);
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int       = ft_sourceinterpolate(cfg, S1bl, template_mri);

cfg    = [];
parcel = ft_sourceparcellate(cfg, source_int, atlas);

% We create a dummy structure where we identify the power values per voxel 
% and use this for subsequent plotting
dummy = atlas;
for i=1:length(parcel.pow)
      dummy.tissue(find(dummy.tissue==i)) = parcel.pow(i);
end

source_int.parcel   = dummy.tissue;
source_int.coordsys = 'mni';
cfg                 = [];
cfg.method          = 'ortho';
cfg.funparameter    = 'parcel';
cfg.funcolormap     = 'jet';
cfg.renderer        = 'zbuffer';
cfg.location        = [-42 -20 6];
cfg.atlas           = atlas;
%cfg.funcolorlim     = [-30 30];
ft_sourceplot(cfg,source_int);

% alternatively
cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'parcel';
%cfg.funcolorlim  = [-30 30];
cfg.funcolormap  = 'jet';
cfg.projmethod   = 'nearest';
cfg.surfinflated = 'surface_inflated_both_caret.mat';
cfg.projthresh   = 0.8;
cfg.camlight     = 'no';
ft_sourceplot(cfg, source_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull


%% without a contrast
%% try new stuff
%--------------
cfg                = [];
cfg.channel        = 'megplanar';
data_preprocessed1 = ft_selectdata(cfg,data_preprocessed);

cfg             = [];
cfg.channel     = 'megplanar';
cfg.grad        = data_preprocessed1.grad;
cfg.headmodel   = headmodel_meg;
cfg.reducerank  = 2;                            % default for MEG is 2, for EEG is 3
cfg.sourcemodel = sourcemodel;
% cfg.normalize   = 'yes';                       % if you are not contrasting
grid            = ft_prepare_leadfield(cfg);

cfg        = [];
cfg.toilim = [-0.5 0];
datapre    = ft_redefinetrial(cfg, data_preprocessed1);
cfg.toilim = [0 0.5];
datapost   = ft_redefinetrial(cfg, data_preprocessed1);

cfg            = [];
cfg.covariance = 'yes';
avgpre         = ft_timelockanalysis(cfg,datapre);
avgpst         = ft_timelockanalysis(cfg,datapost);

figure
imagesc(avgpre.cov)

% reconstruct activity with precomputed filter
cfg                    = [];
cfg.lcmv.projectnoise  = 'yes'; % gives noise estimate
cfg.method             = 'lcmv';
%cfg.method             = 'sam';
cfg.sourcemodel        = grid;
cfg.headmodel          = headmodel_meg;
sourcepre              = ft_sourceanalysis(cfg, avgpre);
sourcepst              = ft_sourceanalysis(cfg, avgpst);

sourcepre.avg.pow = sourcepre.avg.pow ./ sourcepre.avg.noise;
sourcepst.avg.pow = sourcepst.avg.pow ./ sourcepst.avg.noise;

templatefile     = 'Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\template\anatomy\single_subj_T1.nii';
template_mri     = ft_read_mri(templatefile);
sourcepre.pos    = template_grid.pos; % template grid position are in mni space
sourcepst.pos    = template_grid.pos;
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
sourcepre_int    = ft_sourceinterpolate(cfg, sourcepre, template_mri);
sourcepst_int    = ft_sourceinterpolate(cfg, sourcepst, template_mri);

% look at plot
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.location      = [64 -32 8];
cfg.funcolormap   = 'jet';
ft_sourceplot(cfg,sourcepre_int);
ft_sourceplot(cfg,sourcepst_int);

%% other new stuff
templatefile     = 'Z:/Software/MEG_EEG_Toolboxen/fieldtrip-20191127/external/spm8/templates/T1.nii';
template_mri     = ft_read_mri(templatefile);
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int       = ft_sourceinterpolate(cfg, S1bl, template_mri);

cfg    = [];
parcel = ft_sourceparcellate(cfg, source_int, atlas);

% We create a dummy structure where we identify the power values per voxel 
% and use this for subsequent plotting
dummy = atlas;
for i=1:length(parcel.pow)
      dummy.tissue(find(dummy.tissue==i)) = parcel.pow(i);
end

source_int.parcel   = dummy.tissue;
source_int.coordsys = 'mni';
cfg                 = [];
cfg.method          = 'ortho';
cfg.funparameter    = 'parcel';
cfg.funcolormap     = 'jet';
cfg.renderer        = 'zbuffer';
cfg.location        = [-42 -20 6];
cfg.atlas           = atlas;
%cfg.funcolorlim     = [-30 30];
ft_sourceplot(cfg,source_int);

% alternatively
cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'parcel';
%cfg.funcolorlim  = [-30 30];
cfg.funcolormap  = 'jet';
cfg.projmethod   = 'nearest';
cfg.surfinflated = 'surface_inflated_both_caret.mat';
cfg.projthresh   = 0.8;
cfg.camlight     = 'no';
ft_sourceplot(cfg, source_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull










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
cfg.dataset              = [subjectdata.rawdatadir filesep 'empty_pre_sss.fif'];
cfg.trialfun             = 'ft_trialfun_general'; 
cfg.trialdef.triallength = 10;
cfg.trialdef.ntrials     = inf;
cfg                      = ft_definetrial(cfg);
cfg.channel              = 'meg';

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











end % subjects

disp('do_beamforming finished')

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,data_preprocessed);
suptitle('average magnetometer')

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
%cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,data_preprocessed);
suptitle('average gradiometer')

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
%cfg.layout     = 'neuromag306planar.lay';
%cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,data_preprocessed);
suptitle('average gradiometer')

avgpre1 = avgpre;
avgpre1.time = avgpst.time;
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
%cfg.layout     = 'neuromag306mag.lay';
cfg.layout     = 'neuromag306planar.lay';
%cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,avgpre1,avgpst);
suptitle('average gradiometer')

% have a look
% before ica
% cfg                              = [];
% cfg.channel                      = 'megplanar';
% cfg.dataset                      = [subjectdata.rawdatadir filesep filenames{1} '.fif'];
% cfg.blocksize                    = 100; % pulses are in the first 100sec
% cfg.artfctdef.threshold.artifact = artifact_threshold1; % use thresholds for correct file
% cfg.plotevents                   = 'no';
% ft_databrowser(cfg)

% after ica
% cfg                         = [];
% cfg.channel                 = 'megplanar';
% cfg.continuous              = 'yes';
% cfg.blocksize               = 100;
% cfg.artfctdef.threshold.artifact = artifact_threshold1;
% cfg.plotevents              = 'no';
% ft_databrowser(cfg,data_preprocessed) % load preprocessed data with applied ica

