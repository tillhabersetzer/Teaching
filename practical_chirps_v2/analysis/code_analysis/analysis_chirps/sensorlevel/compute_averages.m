close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
% subjects = {'subject04'};
for i = 3:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjects{i-2} = [subject,num2str(i)]; 
end

% check results in plots
check = 0;

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [1, 40];

% apply ica
ica_status = 1;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% apply baseline correction
baseline_correction_status = 1;

% downsample data
downsample_status = 1;
fs_down           = 200;

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

cfg = [];
cfg.channel          = 'meg';
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg                  = ft_timelockanalysis(cfg, data_preprocessed);

%% Noise-covariance estimation
% for a correct noise-covariance estimation it is important that 
% you used the cfg.demean = 'yes';

filenames_noise = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
noise           = give_noise(filenames_noise,subjectdata);

% if config.ica_status
%    noise = reject_independent_components(noise,subjectdata,config.ica_files); 
% end

cfg                  = [];
cfg.channel          = 'meg';
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg,noise);

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
kappa_noise = give_kappa_value(avg_noise.cov,avg_noise.label,{'megmag','megplanar'});
kappa_data  = give_kappa_value(avg.cov,avg.label,{'megmag','megplanar'});

kappa_data(kappa_data<50) = [];
kappa                     = min([kappa_noise,kappa_data]); 
if kappa < 50; error('unexpected kappa value'); end

cfg                 = [];
cfg.channel         = 'meg';
cfg.kappa           = kappa; % ensures use of regularized inverse

if ica_status % necessary for ft_denoise to work
    [i1,i2] = match_str(data_preprocessed.grad.label,avg_noise.grad.label);
    data_preprocessed.grad.chanunit(i1) = avg_noise.grad.chanunit(i2);
    data_preprocessed.grad.chantype(i1) = avg_noise.grad.chantype(i2);
end
data_preprocessed_w = ft_denoise_prewhiten(cfg, data_preprocessed, avg_noise);

cfg      = [];
avg_w    = ft_timelockanalysis(cfg, data_preprocessed_w);

% check averaged signals
%-----------------------
if check
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg,avg_w);
    sgtitle([subjects{s},': all channels'])    
end

if save_data
% save data
%----------
info                  = config;
info.number_of_epochs = length(data_preprocessed.trial);
info.kappa_noise      = kappa_noise;
info.kappa_data       = kappa_data;
info.kappa            = kappa;
info.kappa_order      = {'megmag','megplanar'};

save(fullfile(subjectdata.chirps_sensorlevel,[subjectdata.subjectname,'_averages',add,'.mat']),...
     'info','avg_w','avg','data_preprocessed');
end



end

%% Clean up
addpath(fullfile('Z:','analysis','subject_files'));
addpath(fullfile('Z:','analysis','preprocessing_batch','helper_functions'));
addpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));