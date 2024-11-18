close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about beamforming results
% on a volume grid model.
% It is constructed to contrast the data with the baseline period.
%
% If you are not interested in a contrast, you can normalize the the
% leadfield.
%
% with prewhitening!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
% subjects = {'subject04'};
for i = 3:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjects{i-2} = [subject,num2str(i)]; 
end

% channels of interest for sourcemodelling
coi = 'meg'; % 'megmag','megplanar'

% choose files
% files2preproc = 'stories_maxfilter';

% bandpass fequency
% bp_freq = [4, 30];

% apply ica
ica_status = 1;

% choose files for detected ica components
% ica_files = 'all_maxfilter';

% apply baseline correction
% baseline_correction_status = 0;

% downsample data
% downsample_status = 0;
% fs_down           = 200;

% save data 
save_data = 1;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'));
% addpath for ica functions
addpath(fullfile('Z:','analysis','preprocessing_batch','helper_functions'));
% addpath for preprocessing function
addpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));

if ica_status
    add = '_ica';
else 
    add = '';
end

%% analysis

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj

%% preprocessing
% config.subject                    = subjects{s};
% config.files2preproc              = files2preproc;
% config.bp_freq                    = bp_freq;
% config.ica_status                 = ica_status;
% config.ica_files                  = ica_files;
% config.baseline_correction_status = baseline_correction_status;
% config.downsample_status          = downsample_status;
% config.fs_down                    = fs_down;
% 
% data_preprocessed = preprocess_chirps(config);

eval(subjects{s})
data              = importdata(fullfile(subjectdata.chirps_sensorlevel,[subjectdata.subjectname,'_averages',add,'.mat']));
data_preprocessed = data.data_preprocessed;
config            = data.info;
clear data;

%% noise-covariance estimation
% for a correct noise-covariance estimation it is important that 
% you used the cfg.demean = 'yes';

filenames_noise = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
noise           = give_noise(filenames_noise,subjectdata);

% if config.ica_status
%    noise = reject_independent_components(noise,subjectdata,config.ica_files); 
% end

cfg                  = [];
cfg.channel          = coi;
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg,noise);
avg                  = ft_timelockanalysis(cfg,data_preprocessed);

%% prewhitening
% the following lines detect the location of the first large 'cliff' in the 
% singular value spectrum of the grads and mags
switch coi
    case 'meg'
        sensors = {'megmag','megplanar'};
    case {'megplanar','megmag'}
        sensors = {coi};
end
kappa_noise = give_kappa_value(avg_noise.cov,avg_noise.label,sensors);
kappa_data  = give_kappa_value(avg.cov,avg.label,sensors);

kappa_data(kappa_data<50) = [];
kappa                     = min([kappa_noise,kappa_data]); 
if kappa < 50; error('unexpected kappa value'); end

if ica_status % necessary for ft_denoise to work
    [i1,i2] = match_str(data_preprocessed.grad.label,avg_noise.grad.label);
    data_preprocessed.grad.chanunit(i1) = avg_noise.grad.chanunit(i2);
    data_preprocessed.grad.chantype(i1) = avg_noise.grad.chantype(i2);
end

cfg                 = [];
cfg.channel         = coi;
cfg.kappa           = kappa; % ensures use of regularized inverse
data_preprocessed_w = ft_denoise_prewhiten(cfg, data_preprocessed, avg_noise);
noise_w             = ft_denoise_prewhiten(cfg, noise, avg_noise);

% compute covariance over the complete teamwindow for beamformer
cfg                  = [];
cfg.removemean       = 'yes'; % default for covariance computation
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_w                = ft_timelockanalysis(cfg,data_preprocessed_w);

%% Condition contrast
cfg        = [];
cfg.toilim = [-0.45 0]; % does not make problems with sample size e.g. 500 vs 501
datapre_w  = ft_redefinetrial(cfg, data_preprocessed_w);
cfg.toilim = [0 0.45];
datapst_w  = ft_redefinetrial(cfg, data_preprocessed_w);

cfg                  = [];
cfg.channel          = coi;
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avgpre_w             = ft_timelockanalysis(cfg,datapre_w);
avgpst_w             = ft_timelockanalysis(cfg,datapst_w);

%% downsample data
% after covariance estimation
% cov_all        = avg_w.cov;
% cov_pre        = avgpre_w.cov;
% cov_pst        = avgpst_w.cov;
% 
% cfg            = [];
% cfg.resamplefs = fs_down;
% cfg.detrend    = 'no';
% avg_w          = ft_resampledata(cfg,avg_w);
% avgpre_w       = ft_resampledata(cfg,avgpre_w);
% avgpst_w       = ft_resampledata(cfg,avgpst_w);
% 
% avg_w.cov    = cov_all;
% avgpre_w.cov = cov_pre;
% avgpst_w.cov = cov_pst;     

%% Forward solution
% load headmodel
headmodel        = importdata(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_headmodel.mat'])); % mm
% load sourcemodel 
sourcemodel_grid = importdata(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_grid.mat'])); % mm
sourcemodel_surf = importdata(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_4k.mat'])); % mm
% prepare leadfield
%------------------
cfg                       = [];
cfg.grad                  = data_preprocessed_w.grad;  % sensor information
cfg.channel               = coi;                       % the used channels
cfg.sourcemodel           = sourcemodel_grid;          % source points
cfg.headmodel             = headmodel;                 % volume conduction model
cfg.singleshell.batchsize = 5000;                      % speeds up the computation
%cfg.normalize            = 'yes';                     % if you are not contrasting
leadfield_grid            = ft_prepare_leadfield(cfg); % NOTE: input of the whitened data ensures the correct sensor definition to be used.
cfg.sourcemodel           = sourcemodel_surf; 
leadfield_surf            = ft_prepare_leadfield(cfg);

%% Inverse Solution
kappa_data                = give_kappa_value(avg_w.cov,avg_w.label,sensors);
kappa_data(kappa_data<50) = [];
kappa                     = min([kappa_noise,kappa_data]); 
if kappa < 50; error('unexpected kappa value'); end

% loop over fixed/freely oriented dipoles
%----------------------------------------
fixed_oriented = {'yes','no'};
label          = {'fixed_oriented','freely_oriented'};
for o = 1 % 
cfg                    = [];
cfg.method             = 'lcmv';
%cfg.lcmv.projectnoise = 'yes'; % gives noise estimate for power?
cfg.lcmv.kappa         = kappa;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.fixedori      = fixed_oriented{o}; % project the leadfield onto the orientation of maximum power
cfg.lcmv.weightnorm    = 'unitnoisegain';      % weight normalization prevents the 'center of the head' artifact -> NECESSARY!
cfg.headmodel          = headmodel;
cfg.sourcemodel        = leadfield_grid;
source_grid            = ft_sourceanalysis(cfg,avg_w);
cfg.sourcemodel.filter       = source_grid.avg.filter;
cfg.sourcemodel.filterdimord = source_grid.avg.filterdimord;
sourcepre_grid               = ft_sourceanalysis(cfg,avgpre_w);
sourcepst_grid               = ft_sourceanalysis(cfg,avgpst_w);

cfg.sourcemodel        = leadfield_surf;
source_surf            = ft_sourceanalysis(cfg,avg_w);
cfg.sourcemodel.filter       = source_surf.avg.filter;
cfg.sourcemodel.filterdimord = source_surf.avg.filterdimord;
sourcepre_surf               = ft_sourceanalysis(cfg,avgpre_w);
sourcepst_surf               = ft_sourceanalysis(cfg,avgpst_w);

% save data
%----------
    if save_data
        info                  = config;
        info.prewhitening     = 'yes';
        info.number_of_epochs = length(data_preprocessed.trial);
        info.kappa_noise      = kappa_noise;
        info.kappa_data       = kappa_data;
        info.kappa            = kappa;
        info.kappa_order      = sensors;

        sources_grid.source    = source_grid;
        sources_grid.sourcepre = sourcepre_grid;
        sources_grid.sourcepst = sourcepst_grid;
        sources_grid.info      = info;
        sources_surf.source    = source_surf;
        sources_surf.sourcepre = sourcepre_surf;
        sources_surf.sourcepst = sourcepst_surf;
        sources_surf.info      = info;

%         save(fullfile(subjectdata.chirps_beamformer,[subjectdata.subjectname,'_beamformer_gridmodel_',label{o},'.mat']),...
%         'sources_grid');
%         save(fullfile(subjectdata.chirps_beamformer,[subjectdata.subjectname,'_beamformer_surfacemodel_',label{o},'.mat']),...
%         'sources_surf');

        save(fullfile(subjectdata.chirps_beamformer,[subjectdata.subjectname,'_beamformer_gridmodel',add,'.mat']),...
        'sources_grid');
        save(fullfile(subjectdata.chirps_beamformer,[subjectdata.subjectname,'_beamformer_surfacemodel',add,'.mat']),...
        'sources_surf');

    end
    clear source_grid sourcepre_grid sourcepst_grid source_surf sourcepre_surf sourcespst_suf ...
          info sources_surf sources_grid
end

end

%% Clean up
rmpath(fullfile('Z:','analysis','subject_files'));
rmpath(fullfile('Z:','analysis','preprocessing_batch','helper_functions'));
rmpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));