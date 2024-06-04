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
subjects = {'subject04'};

% channels of interest for sourcemodelling
coi = 'meg';

% check results in plots
check = 0;

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [4, 30];

% apply ica
ica_status = 0;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% apply baseline correction
baseline_correction_status = 1;

% downsample data
downsample_status = 0;
fs_down           = 300;
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
% change paths if you work with Linux
subjectdata.rawdatadir = replace(subjectdata.rawdatadir,old,new);    
    
%% Averaging and noise-covariance estimation
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

cfg                  = [];
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = [-inf 0]; % timepoints before the zero timepoints in trials
avg                  = ft_timelockanalysis(cfg, data_preprocessed);

% noise covariance comparison
%----------------------------
if check
    figure
    subplot(2,2,1)
    imagesc(avg_noise.cov);
    title('empty room measurements')
    subplot(2,2,2)
    imagesc(avg.cov);
    title('baseline interval')
    sgtitle([subjectdata.subjectname,': noise covariance matrix'])
    
    selmag  = ft_chantype(avg_noise.label, 'megmag');
    selgrad = ft_chantype(avg_noise.label, 'megplanar');
    C = avg_noise.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    subplot(2,2,3)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('empty room measurements')
    
    selmag  = ft_chantype(avg.label, 'megmag');
    selgrad = ft_chantype(avg.label, 'megplanar');
    C = avg.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    subplot(2,2,4)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);  
    title('baseline interval')
    
    % singular values
    %----------------
    [~,sv,~] = svd(avg_noise.cov);
    figure
    subplot(1,2,1)
    plot(log10(diag(sv)),'o');
    title('empty room measurements')
    [~,sv,~] = svd(avg.cov);
    subplot(1,2,2)
    plot(log10(diag(sv)),'o');
    title('baseline interval')
end

%% Forward solution
% load headmodel
headmodel = importdata(replace(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_headmodel.mat']),old,new)); % mm
% load sourcemodel (4k is already enough)
sourcemodel = importdata(replace(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_4k.mat']),old,new)); % mm

headmodel   = ft_convert_units(headmodel, data_preprocessed.grad.unit);
sourcemodel = ft_convert_units(sourcemodel, data_preprocessed.grad.unit);

% prepare leadfield
%------------------
cfg                       = [];
cfg.grad                  = data_preprocessed.grad;    % sensor information
cfg.channel               = coi;                       % the used channels
cfg.sourcemodel           = sourcemodel;               % source points
cfg.headmodel             = headmodel;                 % volume conduction model
cfg.singleshell.batchsize = 5000;                      % speeds up the computation
leadfield                 = ft_prepare_leadfield(cfg);

% check coordinate systems
%-------------------------
if check
    dataset = replace([subjectdata.rawdatadir filesep 'fixer_olsa_tsss_mc.fif'],old,new);
    shape   = ft_read_headshape(dataset,'unit',data_preprocessed.grad.unit); 

    figure
    ft_plot_headmodel(headmodel,'facealpha',0.5);
    hold on
    ft_plot_sens(data_preprocessed.grad, 'style', '*b');
    ft_plot_headshape(shape);
    ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', ...
                 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);
end
%% Inverse solution
% modify noise covariance estimation
avg.cov = avg_noise.cov;

cfg                    = [];
cfg.method             = 'mne';
cfg.sourcemodel        = leadfield;
cfg.headmodel          = headmodel;
cfg.mne.prewhiten      = 'yes'; % combine different sensortypes with different units (megmag,megplanar)
cfg.mne.lambda         = 3;     % scaling factor for noise
cfg.mne.scalesourcecov = 'yes';
cfg.mne.keepfilter     = 'yes';
source                 = ft_sourceanalysis(cfg,avg);

% save data
%----------
save(replace(fullfile(subjectdata.chirps_mne,[subjectdata.subjectname,'_mne_sources.mat']),old,new),'source')

info.number_of_epochs = length(data_preprocessed.trial);
save(replace(fullfile(subjectdata.chirps_mne,[subjectdata.subjectname,'_mne_sources_info.mat']),old,new),'info')

% fid = fopen(replace(fullfile(subjectdata.chirps_mne,[subjectdata.subjectname,'mne_sources_info.txt']),old,new),'wt');
% fprintf(fid, 'number of epochs left: %s',num2str(length(data_preprocessed.trial)));
% fclose(fid);

%source.inside = mask;
%source.coordsys = 'neuromag';
%% Visualize
if check 
    signal = source.avg.pow(:,640); % 140ms
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
    %ft_sourcemovie(cfg,sd);
    ft_sourcemovie(cfg,source); % without projection -> looks better?
end

end

% use atlas with mask
% or this labels (more accurate) '.R.aparc.a2009s.4k_fs_LR.label.gii'
% path_label_r = replace(fullfile(subjectdata.sourcemodel,'workbench',[subjectdata.subjectname,'.R.aparc.4k_fs_LR.label.gii']),old,new);
% path_mesh_r  = replace(fullfile(subjectdata.sourcemodel,'workbench',[subjectdata.subjectname,'.R.midthickness.4k_fs_LR.surf.gii']),old,new);
% 
% atlas_r = ft_read_atlas({path_label_r,path_mesh_r});
% atlas_l = ft_read_atlas({strrep(path_label_r,'.R.','.L.'),strrep(path_mesh_r,'.R.','.L.')});
% 
% roi  = {'L_superiortemporal','R_superiortemporal'};
% mask = mne_generate_mask(roi,atlas_l,atlas_r,sourcemodel);

%% Clean up
rmpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new))
rmpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new))
rmpath(replace(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions'],old,new));
