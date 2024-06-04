 close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the averages 
% - for prewhitened data and not-prewhitened (raw) data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject04'};

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
downsample_status = 1;
fs_down           = 200;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);
% addpath for ica functions
addpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);
% addpath for preprocessing function
addpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);

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
avg = ft_timelockanalysis(cfg, data_preprocessed);

cfg                = [];
cfg.method         = 'sum';
cfg.updatesens     = 'yes'; 
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.5 0];
avg_cmb            = ft_combineplanar(cfg,avg);
    
%% Noise-covariance estimation
% for a correct noise-covariance estimation it is important that 
% you used the cfg.demean = 'yes';

% decide which noise to use
%--------------------------
filenames_noise = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
noise           = give_noise(filenames_noise,subjectdata);

% or baseline
%------------
% cfg        = [];
% cfg.toilim = [-0.5,0];
% noise      = ft_redefinetrial(cfg, data_preprocessed);

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
    [~,vs,~] = svd(avg_noise.cov);
    subplot(1,2,2)
    plot(log10(diag(vs)),'o');
    grid on
    title('Singular values of a MEG sensor covariance matrix')
    sgtitle(subjectdata.subjectname)
end

%% prewhitening
% the following lines detect the location of the first large 'cliff' in the 
% singular value spectrum of the grads and mags
kappa_noise = give_kappa_value(avg_noise.cov,avg_noise.label,{'megmag','megplanar'});
kappa_data  = give_kappa_value(avg.cov,avg.label,{'megmag','megplanar'});

cfg                 = [];
cfg.channel         = 'meg';
cfg.kappa           = min([kappa_noise,kappa_data]); % ensures use of regularized inverse
data_preprocessed_w = ft_denoise_prewhiten(cfg, data_preprocessed, avg_noise);

cfg            = [];
cfg.covariance = 'yes';
avg_w          = ft_timelockanalysis(cfg, data_preprocessed_w);

cfg                = [];
cfg.method         = 'sum';
cfg.updatesens     = 'yes'; 
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.5 0];
avg_cmb_w          = ft_combineplanar(cfg,avg_w);

% check if whitening and average commutate - are the same!
%---------------------------------------------------------
cfg              = [];
cfg.channel      = 'meg';
cfg.kappa        = min(kappa_noise); % ensures use of regularized inverse
avg_w_afterwards = ft_denoise_prewhiten(cfg, avg, avg_noise);

disp(isequal(avg_w.avg,avg_w_afterwards.avg));
disp(max(abs(avg_w.avg)-abs(avg_w_afterwards.avg),[],'all'));

if check
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg,avg_w,avg_w_afterwards);
    sgtitle([subjects{s},': all channels'])   
end

% Check Global Field Power
%-------------------------
selmag   = ft_chantype(avg_w.label, 'megmag');
selgrad  = ft_chantype(avg_w.label, 'megplanar');
gfp_mag  = sum(avg_w.avg(selmag,:).^2,1)./kappa_noise(1);
gfp_grad = sum(avg_w.avg(selgrad,:).^2,1)./kappa_noise(2);

if check
    figure
    subplot(1,2,1)
    plot(avg_w.time,gfp_mag)
    title('magnetometer')
    subplot(1,2,2)
    plot(avg_w.time,gfp_grad)
    title('gradiometer')
    sgtitle('global field power')
end

%% Plot some data - timecourses

% magnetometer
%-------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
ft_multiplotER(cfg,avg);
suptitle('averaged magnetometer')

figure
ft_multiplotER(cfg,avg_w);
suptitle('prewhitened average magnetometer')

% gradiometer
%-------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
ft_multiplotER(cfg,avg);
suptitle('averaged gradiometer')

figure
ft_multiplotER(cfg,avg_w);
suptitle('prewhitened average gradiometer')

% combined gradiometer
%---------------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg,avg_cmb);
suptitle('averaged gradiometer combined')

figure
ft_multiplotER(cfg,avg_cmb_w);
suptitle('prewhitened averaged gradiometer combined')

%% Plot some data - topoplots

% magnetometer
%-------------
figure
cfg          = [];
cfg.xlim     = [0.07, 0.140];
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306mag.lay';
ft_topoplotER(cfg,avg);
suptitle('averaged magnetometer')

figure
ft_topoplotER(cfg,avg_w);
suptitle('prewhitened average magnetometer')

% gradiometer
%-------------
figure
cfg          = [];
cfg.xlim     = [0.07, 0.140];
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306planar.lay';
ft_topoplotER(cfg,avg);
suptitle('averaged gradiometer')

figure
ft_topoplotER(cfg,avg_w);
suptitle('prewhitened average gradiometer')

% combined gradiometer
%---------------------
figure
cfg          = [];
cfg.xlim     = [0.07, 0.140];
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306cmb.lay';
ft_topoplotER(cfg,avg_cmb);
suptitle('averaged gradiometer combined')

figure
ft_topoplotER(cfg,avg_cmb_w);
suptitle('prewhitened averaged gradiometer combined')

%% Plot some data - topoplot timeseries

% magnetometer
%-------------
figure
cfg          = [];
cfg.xlim     = [-0.1 : 0.1 : 0.5];  
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306mag.lay';
ft_topoplotER(cfg,avg);
suptitle('averaged magnetometer')

figure
ft_topoplotER(cfg,avg_w);
suptitle('prewhitened average magnetometer')

% gradiometer
%-------------
figure
cfg          = [];
cfg.xlim     = [-0.1 : 0.1 : 0.5];  
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306planar.lay';
ft_topoplotER(cfg,avg);
suptitle('averaged gradiometer')

figure
ft_topoplotER(cfg,avg_w);
suptitle('prewhitened average gradiometer') 

% combined gradiometer
%---------------------
figure
cfg          = [];
cfg.xlim     = [-0.1 : 0.1 : 0.5];  
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306cmb.lay';
ft_topoplotER(cfg,avg_cmb);
suptitle('averaged gradiometer combined')

figure
ft_topoplotER(cfg,avg_cmb_w);
suptitle('prewhitened averaged gradiometer combined')

end % subjects

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files']);
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);
rmpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);