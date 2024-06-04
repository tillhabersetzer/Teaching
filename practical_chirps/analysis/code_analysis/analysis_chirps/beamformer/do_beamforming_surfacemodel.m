 close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about beamforming results
% on a surface based description. 
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
%subjects = {'subject03','subject04','subject05','subject06','subject07'};
subjects = {'subject04'};

% channels of interest for sourcemodelling
coi = 'meg'; % 'megmag','megplanar'

% check results in plots
check = 0;

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [1, 40];

% apply ica
ica_status = 0;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% apply baseline correction
baseline_correction_status = 1;

% downsample data
downsample_status = 1;
fs_down           = 200;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'));
% addpath for ica functions
addpath(fullfile('Z:','analysis','preprocessing_batch','helper_functions'));
% addpath for preprocessing function
addpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));

%% analysis

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj

%% preprocessing
config.subject                    = subjects{s};
config.files2preproc              = files2preproc;
config.bp_freq                    = bp_freq;
config.ica_status                 = ica_status;
config.ica_files                  = ica_files;
config.baseline_correction_status = baseline_correction_status;
config.downsample_status          = downsample_status;
config.fs_down                    = fs_down;

data_preprocessed = preprocess_chirps(config);

eval(subjects{s})

%% noise-covariance estimation
% for a correct noise-covariance estimation it is important that 
% you used the cfg.demean = 'yes';

filenames_noise = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
noise           = give_noise(filenames_noise,subjectdata);

cfg                  = [];
cfg.channel          = coi;
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg,noise);
avg                  = ft_timelockanalysis(cfg,data_preprocessed);

if check
    % noise covariance matrix
    %------------------------
    selmag  = ft_chantype(avg_noise.label, 'megmag');
    selgrad = ft_chantype(avg_noise.label, 'megplanar');
    C = avg_noise.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    figure
    subplot(1,2,1)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('MEG sensor covariance matrix')
    
    % singular values
    %----------------
    [~,sv,~] = svd(avg_noise.cov);
    subplot(1,2,2)
    plot(log10(diag(sv)),'o');
    grid on
    title('Singular values of a MEG sensor covariance matrix')
    sgtitle(subjectdata.subjectname)
end

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

cfg                 = [];
cfg.channel         = coi;
cfg.kappa           = min([kappa_noise,kappa_data]); % ensures use of regularized inverse
data_preprocessed_w = ft_denoise_prewhiten(cfg, data_preprocessed, avg_noise);
noise_w             = ft_denoise_prewhiten(cfg, noise, avg_noise);

% compute covariance over the complete teamwindow for beamformer
cfg                  = [];
cfg.removemean       = 'yes'; % default for covariance computation
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_w                = ft_timelockanalysis(cfg,data_preprocessed_w);

cfg        = [];
cfg.toilim = [-0.45 0]; % does not make problems with sample size e.g. 500 vs 501
datapre_w  = ft_redefinetrial(cfg, data_preprocessed_w);
cfg.toilim = [0 0.45];
datapst_w = ft_redefinetrial(cfg, data_preprocessed_w);

cfg                = [];
cfg.channel        = coi;
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avgpre_w           = ft_timelockanalysis(cfg,datapre_w);
avgpst_w           = ft_timelockanalysis(cfg,datapst_w);

if check 
    cfg                  = [];
    cfg.channel          = coi;
    cfg.covariance       = 'yes';
    cfg.covariancewindow = 'all';
    avg_noise_w          = ft_timelockanalysis(cfg,noise_w);

    % noise and data covariance matrix after whitening
    %-------------------------------------------------
    selmag  = ft_chantype(avg_noise_w.label, 'megmag');
    selgrad = ft_chantype(avg_noise_w.label, 'megplanar');
    C = avg_noise_w.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    figure
    subplot(2,2,1)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('noise')
    
    selmag  = ft_chantype(avg_w.label, 'megmag');
    selgrad = ft_chantype(avg_w.label, 'megplanar');
    C = avg_w.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    subplot(2,2,2)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('data')
    
    % singular values after whitening
    %--------------------------------
    [~,sv,~] = svd(avg_noise_w.cov);
    subplot(2,2,3)
    plot(log10(diag(sv)),'o');
    grid on
    title('noise')
    
    [~,sv,~] = svd(avg_w.cov);
    subplot(2,2,4)
    plot(log10(diag(sv)),'o');
    grid on
    title('data')
    sgtitle([subjectdata.subjectname,': after whitening'])
end

%% Forward solution
% load headmodel
headmodel = importdata(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_headmodel.mat'])); % mm
% load sourcemodel (4k is already enough)
sourcemodel = importdata(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_4k.mat'])); % mm
% load inflated model for visualization
inflated = importdata(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_inflated_4k.mat'])); % mm

headmodel   = ft_convert_units(headmodel, avg_w.grad.unit);
sourcemodel = ft_convert_units(sourcemodel, avg_w.grad.unit);

% prepare leadfield
%------------------
cfg                       = [];
cfg.grad                  = data_preprocessed_w.grad;  % sensor information
cfg.channel               = coi;                       % the used channels
cfg.sourcemodel           = sourcemodel;               % source points
cfg.headmodel             = headmodel;                 % volume conduction model
cfg.singleshell.batchsize = 5000;                      % speeds up the computation
%cfg.normalize             = 'yes';                     % if you are not contrasting
leadfield                 = ft_prepare_leadfield(cfg); % NOTE: input of the whitened data ensures the correct sensor definition to be used.

% check coordinate systems
%-------------------------
if check
    dataset = fullfile(subjectdata.rawdatadir,'fixer_olsa_tsss_mc.fif');
    shape   = ft_read_headshape(dataset,'unit',avg_w.grad.unit); 

    figure
    ft_plot_headmodel(headmodel,'facealpha',0.5);
    hold on
    ft_plot_sens(ft_convert_units(data_preprocessed.grad,avg_w.grad.unit), 'style', '*b');
    ft_plot_headshape(shape);
    ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', ...
                 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);
end

%% downsample data
% after covariance estimation
cov_all        = avg_w.cov;
cov_pre        = avgpre_w.cov;
cov_pst        = avgpst_w.cov;

cfg            = [];
cfg.resamplefs = fs_down;
cfg.detrend    = 'no';
avg_w          = ft_resampledata(cfg,avg_w);
avgpre_w       = ft_resampledata(cfg,avgpre_w);
avgpst_w       = ft_resampledata(cfg,avgpst_w);

avg_w.cov    = cov_all;
avgpre_w.cov = cov_pre;
avgpst_w.cov = cov_pst; 

%% Inverse Solution
kappa_data            = give_kappa_value(avg_w.cov,avg_w.label,sensors);

cfg                   = [];
cfg.method            = 'lcmv';
%cfg.lcmv.projectnoise = 'yes'; % gives noise estimate for power?
cfg.lcmv.kappa        = min([kappa_noise,kappa_data]);
cfg.lcmv.keepfilter   = 'yes';
cfg.lcmv.fixedori     = 'yes';  % project the leadfield onto the orientation of maximum power
cfg.lcmv.weightnorm   = 'unitnoisegain';
cfg.headmodel         = headmodel;
cfg.sourcemodel       = leadfield;
source                = ft_sourceanalysis(cfg,avg_w);

% Sam does not work?
% cfg                   = [];
% cfg.method            = 'sam';
% %cfg.lcmv.projectnoise = 'yes'; % gives noise estimate for power?
% cfg.sam.kappa        = min(kappa_noise);
% cfg.sam.keepfilter   = 'yes';
% %cfg.sam.fixedori     = 'yes';  % does not exist for sam
% cfg.sam.weightnorm   = 'unitnoisegain';
% cfg.headmodel         = headmodel;
% cfg.sourcemodel       = leadfield;
% source                = ft_sourceanalysis(cfg,avg_w);

if check
    % absolute dipole moments
    %------------------------
    cfg           = [];
    cfg.operation = 'abs';
    cfg.parameter = 'mom';
    sourceabsmom  = ft_math(cfg, source);

    cfg           = [];
    %cfg.clim      = [-2.5 2.5];
    cfg.colormap  = 'parula';
    cfg.parameter = 'mom';
    % replace the original dipole positions with those of the inflated surface
    source.pos    = inflated.pos;
    ft_sourceplot_interactive(cfg, source);
    sourceabsmom.pos = inflated.pos;
    ft_sourceplot_interactive(cfg, sourceabsmom);   
end

%% Condition contrast
cfg                          = [];
cfg.method                   = 'lcmv';
cfg.lcmv.kappa               = min([kappa_noise,kappa_data]);
cfg.lcmv.keepfilter          = 'yes';
cfg.lcmv.fixedori            = 'yes';
cfg.lcmv.weightnorm          = 'unitnoisegain';
cfg.headmodel                = headmodel;
cfg.sourcemodel              = leadfield;
cfg.sourcemodel.filter       = source.avg.filter;
cfg.sourcemodel.filterdimord = source.avg.filterdimord;
sourcepre                    = ft_sourceanalysis(cfg, avgpre_w);
sourcepst                    = ft_sourceanalysis(cfg, avgpst_w);

% change time for subtracting data
sourcepre.time = sourcepst.time;

cfg             = [];
cfg.operation   = 'abs';
cfg.parameter   = 'mom';
sourcepreabsmom = ft_math(cfg, sourcepre);
sourcepstabsmom = ft_math(cfg, sourcepst);

cfg           = [];
cfg.parameter = 'mom';
%cfg.operation = '((x1-x2)./x2)*100'; %percentage change from baseline
cfg.operation = 'subtract';
source_diff   = ft_math(cfg,sourcepstabsmom,sourcepreabsmom);

if check
    cfg           = [];
    %cfg.clim      = [-2.5 2.5];
    cfg.colormap  = 'parula';
    cfg.parameter = 'mom';
    cfg.has_diff  = true; % only when you used subtract method
    sourcepreabsmom.pos = inflated.pos;
    sourcepstabsmom.pos = inflated.pos;
    source_diff.pos     = inflated.pos;
    ft_sourceplot_interactive(cfg, sourcepreabsmom,sourcepstabsmom,source_diff); 
end

% save data
%----------
info                  = config;
info.prewhitening     = 'yes';
info.number_of_epochs = length(data_preprocessed.trial);
info.kappa_noise      = kappa_noise;
info.kappa_data       = kappa_data;
info.kappa_order      = sensors;

save(fullfile(subjectdata.chirps_beamformer,[subjectdata.subjectname,'_beamformer_surfacemodel.mat']),...
    'source','sourcepre','sourcepst','info');

end

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])
rmpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);