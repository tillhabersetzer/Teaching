close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!SCRIPT WORKS ALSO WITH LINUX!!!!
% to work on windows replace new = '/media/till/Samsung_T5' with new = old
% then everthing with the pays stays the same
%
% This script tries to give you a first overview about the averages and
% minimum norm estimates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
%subjects = {'subject19'};
for i = 3:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjects{i-2} = [subject,num2str(i)]; 
end

% channels of interest for sourcemodelling
coi = 'meg';

% check results in plots
check = 0;

% choose files
% files2preproc = 'stories_maxfilter';

% bandpass fequency
% bp_freq = [1, 40];

% apply ica
ica_status =1;

% choose files for detected ica components
% ica_files = 'all_maxfilter';

% apply baseline correction
% baseline_correction_status = 1;

% downsample data
% downsample_status = 1;
% fs_down           = 200;

% save data
save_data = 1;
%--------------------------------------------------------------------------

% addpath for subject_files information
old = 'Z:'; 
%new = '/media/till/Samsung_T5';
new = 'Z:';

% addpath for subject_files information
addpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new));
% addpath for ica functions
addpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new));
% addpath for preprocessing function
addpath(replace(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions'],old,new));

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
% change paths if you work with Linux
subjectdata.rawdatadir = replace(subjectdata.rawdatadir,old,new);  
data                   = importdata(fullfile(subjectdata.chirps_sensorlevel,[subjectdata.subjectname,'_averages',add,'.mat']));
data_preprocessed      = data.data_preprocessed;
config                 = data.info;
clear data;
    
%% Noise-covariance estimation
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


cfg                  = [];
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise_w          = ft_timelockanalysis(cfg,noise_w);

% for contrasts
%--------------
cfg        = [];
cfg.toilim = [-0.45 0]; % does not make problems with sample size e.g. 500 vs 501
datapre_w  = ft_redefinetrial(cfg, data_preprocessed_w);
cfg.toilim = [0 0.45];
datapst_w  = ft_redefinetrial(cfg, data_preprocessed_w);

cfg         = [];
cfg.channel = coi;
avg_w       = ft_timelockanalysis(cfg, data_preprocessed_w);

avgpre_w   = ft_timelockanalysis(cfg, datapre_w);
avgpst_w   = ft_timelockanalysis(cfg, datapst_w);

% check covariance
%-----------------
if check 
    % noise and data covariance matrix after whitening
    %-------------------------------------------------
    selmag  = ft_chantype(avg_noise_w.label, 'megmag');
    selgrad = ft_chantype(avg_noise_w.label, 'megplanar');
    C = avg_noise_w.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    figure
    subplot(1,2,1)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('noise')
    
    % singular values after whitening
    %--------------------------------
    [~,sv,~] = svd(avg_noise_w.cov);
    subplot(1,2,2)
    plot(log10(diag(sv)),'o');
    grid on
    title('noise')
end

% additional rejection of bad epochs
%-----------------------------------
% cfg        = [];
% cfg.method = 'summary';
% cfg.layout = 'neuromag306mag_helmet.mat';
% dataw_meg  = ft_rejectvisual(cfg, data_preprocessed_w);

%% Forward solution
% load headmodel
headmodel = importdata(replace(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_headmodel.mat']),old,new)); % mm
% load sourcemodel (4k is already enough)
sourcemodel = importdata(replace(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_4k.mat']),old,new)); % mm

headmodel   = ft_convert_units(headmodel, data_preprocessed_w.grad.unit);
sourcemodel = ft_convert_units(sourcemodel, data_preprocessed_w.grad.unit);

% prepare leadfield
%------------------
cfg                       = [];
cfg.grad                  = data_preprocessed_w.grad;  % sensor information
cfg.channel               = coi;                       % the used channels
cfg.sourcemodel           = sourcemodel;               % source points
cfg.headmodel             = headmodel;                 % volume conduction model
cfg.singleshell.batchsize = 5000;                      % speeds up the computation
leadfield                 = ft_prepare_leadfield(cfg); % NOTE: input of the whitened data ensures the correct sensor definition to be used.

% check coordinate systems
%-------------------------
if check
    dataset = replace([subjectdata.rawdatadir filesep 'fixer_olsa_tsss_mc.fif'],old,new);
    shape   = ft_read_headshape(dataset,'unit',data_preprocessed_w.grad.unit); 

    figure
    ft_plot_headmodel(headmodel,'facealpha',0.5);
    hold on
    ft_plot_sens(data_preprocessed.grad, 'style', '*b');
    ft_plot_headshape(shape);
    ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', ...
                 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);
end
%% Inverse solution
% modify noise covariance estimation - necessary for identical filters!
avg_w.cov    = avg_noise_w.cov;
avgpre_w.cov = avg_noise_w.cov;
avgpst_w.cov = avg_noise_w.cov;

cfg                    = [];
cfg.method             = 'mne';
cfg.sourcemodel        = leadfield;
cfg.headmodel          = headmodel;
cfg.mne.prewhiten      = 'yes'; % combine different sensortypes with different units (megmag,megplanar)
cfg.mne.lambda         = 3;     % scaling factor for noise
cfg.mne.scalesourcecov = 'yes';
cfg.mne.keepfilter     = 'yes';
source                 = ft_sourceanalysis(cfg,avg_w);
source_pre             = ft_sourceanalysis(cfg,avgpre_w);
source_pst             = ft_sourceanalysis(cfg,avgpst_w);

if save_data
% save data
%----------
info                  = config;
info.prewhitening     = 'yes';
info.number_of_epochs = length(data_preprocessed.trial);
info.kappa_noise      = kappa_noise;
info.kappa_data       = kappa_data;
info.kappa            = kappa;
info.kappa_order      = sensors;

save(replace(fullfile(subjectdata.chirps_mne,[subjectdata.subjectname,'_mne_sources',add,'.mat']),old,new),...
     'source','source_pre','source_pst','info');
end
% fid = fopen(replace(fullfile(subjectdata.chirps_mne,[subjectdata.subjectname,'mne_sources_info.txt']),old,new),'wt');
% fprintf(fid, 'number of epochs left: %s',num2str(length(data_preprocessed.trial)));
% fclose(fid);

%source.inside = mask;
%source.coordsys = 'neuromag';
%% Visualize
if check 
    %signal = source.avg.pow(:,640); % 140ms
    idx = dsearchn(source.time',0.140);
    signal = source.avg.pow(:,idx); % 140ms
    figure
    ft_plot_mesh(source, 'vertexcolor', signal);
    view([180 0]);
    h = light; 
    set(h, 'position', [0 1 0.2]); 
    lighting gouraud;
    material dull;
    title([subjectdata.subjectname,': 140ms'])

    % movie
    cfg            = [];
    % project the source to its strongest orientation, i.e. the direction that explains most of the source variance. 
    % That is equivalent to taking the largest eigenvector of the source
    % timeseries.
    cfg.projectmom = 'yes';
    sd             = ft_sourcedescriptives(cfg,source);

    cfg              = [];
    cfg.funparameter = 'pow';
    ft_sourcemovie(cfg,source); % without projection -> looks better?
    
    cfg              = [];
    cfg.funparameter = 'pow';
    ft_sourcemovie(cfg,sd);
end

%% parcellation

if check     
% change for linux
subjectdata.sourcemodel = replace(subjectdata.sourcemodel,old,new);  

% generate surface based atlas for both hemispheres
atlasname = '.L.aparc.a2009s.4k_fs_LR.label.gii';
meshname  = '.L.midthickness.4k_fs_LR.surf.gii';
atlas_big = generate_atlas(atlasname,meshname,sourcemodel,subjectdata);
atlas_big = ft_convert_units(atlas_big,data_preprocessed_w.grad.unit);

atlasname = '.L.aparc.4k_fs_LR.label.gii';
meshname  = '.L.midthickness.4k_fs_LR.surf.gii';
atlas_small = generate_atlas(atlasname,meshname,sourcemodel,subjectdata);
atlas_small = ft_convert_units(atlas_small,data_preprocessed_w.grad.unit);

% not necessary - source and atlas point are identical
% cfg              = [];
% cfg.voxelcoord   = 'no';
% cfg.parameter    = 'mom';
% cfg.interpmethod = 'nearest';
% source_int       = ft_sourceinterpolate(cfg,source,atlas);

cfg        = [];
cfg.method = 'mean';
parcel     = ft_sourceparcellate(cfg, source, atlas_big);

cfg              = [];
cfg.funparameter = 'pow';
%ft_sourcemovie(cfg,sd);
ft_sourcemovie(cfg,parcel);
end

end

% roi  = {'L_superiortemporal','R_superiortemporal'};
% mask = mne_generate_mask(roi,atlas_l,atlas_r,sourcemodel);

%% Clean up
rmpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new))
rmpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new))
rmpath(replace(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions'],old,new));
