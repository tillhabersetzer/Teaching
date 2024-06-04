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
baseline_correction_status = 0;

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
    
%% Noise-covariance estimation
% for a correct noise-covariance estimation it is important that 
% you used the cfg.demean = 'yes';

filenames_noise = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
noise           = give_noise(filenames_noise,subjectdata);

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
    [~,s,~] = svd(avg_noise.cov);
    subplot(1,2,2)
    plot(log10(diag(s)),'o');
    grid on
    title('Singular values of a MEG sensor covariance matrix')
    sgtitle(subjectdata.subjectname)
end

%% prewhitening
% the following lines detect the location of the first large 'cliff' in the 
% singular value spectrum of the grads and mags
kappa_noise = give_kappa_value(avg_noise.cov,avg_noise.label,{'megmag','megplanar'});

cfg                 = [];
cfg.channel         = 'meg';
cfg.kappa           = min(kappa_noise); % ensures use of regularized inverse
data_preprocessed_w = ft_denoise_prewhiten(cfg, data_preprocessed, avg_noise);

cfg        = [];
cfg.toilim = [-0.45 0];
datapre_w  = ft_redefinetrial(cfg,data_preprocessed_w);
cfg.toilim = [0 0.45];
datapst_w  = ft_redefinetrial(cfg,data_preprocessed_w);

datapre_w.time = datapst_w.time;

% optional - test against averaged baseline
%------------------------------------------
% average over time (columns) and replace matrix data with average
datapre_w.trial = cellfun(@(x) repmat(mean(x,2),1,size(x,2)),datapre_w.trial,'UniformOutput',false);

cfg      = [];
avgpre_w = ft_timelockanalysis(cfg, datapre_w);
avgpst_w = ft_timelockanalysis(cfg, datapst_w);

cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffect_w   = ft_math(cfg,avgpst_w, avgpre_w);

% check averaged signals
%-----------------------
if check
    % without baseline corrected data not useful
    cfg            = [];
    cfg.method     = 'sum';
    cfg.updatesens = 'yes'; 
    avgpre_cmb_w   = ft_combineplanar(cfg,datapre_w);
    avgpst_cmb_w   = ft_combineplanar(cfg,datapst_w);

    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg,avgpre_w,avgpst_w);
    sgtitle('all channels')    
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb_helmet.mat';
    ft_multiplotER(cfg,avgpre_cmb_w,avgpst_cmb_w);
    sgtitle('combined gradiometers')    
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg, raweffect_w);
    sgtitle('all channels - raweffect')
end

%% statistic

cfg                  = [];
cfg.channel          = 'meg';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 4;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb           = [];
cfg_neighb.channel   = 'meg';
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, datapre_w);
%ft_neighbourplot(cfg,datapre_w)

ntrials                       = length(datapre_w.trial);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable 
cfg.uvar   = 2; % row of design matrix that contains unit of observation

stat = ft_timelockstatistics(cfg, datapst_w, datapre_w);

% make nice plots
%----------------
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust  = find(pos_cluster_pvals < stat.cfg.alpha); % use all significant clusters
pos               = ismember(stat.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < stat.cfg.alpha);
neg               = ismember(stat.negclusterslabelmat, neg_signif_clust);

timestep      = 0.025; 
sampling_rate = datapre_w.fsample; 
sample_count  = length(stat.time);
j = [0:timestep:0.45];                       % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in samples

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect_w.label, stat.label);

if check
% Plot over raweffect
%--------------------
figure
for k = 1:18
   subplot(4,5,k);
   cfg      = [];
   cfg.xlim = [j(k) j(k+1)];  
   pos_int     = zeros(numel(raweffect_w.label),1);
   neg_int     = zeros(numel(raweffect_w.label),1);
   %pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   %neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
   % less severe
   pos_int(i1) = any(pos(i2, m(k):m(k+1)), 2);
   neg_int(i1) = any(neg(i2, m(k):m(k+1)), 2);
   cfg.highlight        = 'on';
   % Get the index of each significant channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 10;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306all_helmet.mat';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg, raweffect_w);
end

% Plot over t-values
%-------------------
% change t-values to abs for plotting
stati      = stat;
stati.stat = abs(stat.stat);

figure
for k = 1:18
   subplot(4,5,k);
   cfg      = [];
   cfg.xlim = [j(k) j(k+1)];  
   %cfg.zlim = [0,5];
   cfg.zlim = [0 max(stati.stat,[],'all')];
   cfg.highlight        = 'on';
   % int = all(stati.mask(:, m(k):m(k+1)), 2);
   % less severe
   int = any(stati.mask(:, m(k):m(k+1)), 2); 
   cfg.highlightchannel = find(int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 10;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306all_helmet.mat';
   cfg.interactive      = 'no';
   cfg.parameter        = 'stat';
   ft_topoplotER(cfg, stati);
   set(gca,'ColorScale','log') % if you want to plot on a logarithmic scale
end
c           = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'t-val');


cfg               = [];
cfg.layout        = 'neuromag306all_helmet.mat';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor    = 'r';
figure;
ft_multiplotER(cfg,stat);

% plot N1 and P2 topography in units of t-value
cfg           = [];
cfg.layout    = 'neuromag306all_helmet.mat';
cfg.parameter = 'stat';
%cfg.zlim      = [0 max(stati.stat,[],'all')];
cfg.zlim      = [0,5];
%cfg.zlim      = [-5,5];
cfg.colormap  = parula;
cfg.marker    = 'off';
cfg.style     = 'fill';
cfg.comment   = 'off';
cfg.colorbar  = 'no';

% time windows for topoplots
windows = [0,0.05;0.05,0.075;0.075,0.1;0.1,0.125;0.125,0.15;0.15,0.175;0.175,0.2;...
           0.2,0.225;0.225,0.25;0.25,0.3];
figure
sgtitle('topography')
for i = 1:10
cfg.xlim = windows(i,:);
subplot(2,5,i)
ft_topoplotER(cfg,stati);
title([num2str(cfg.xlim(1)),'-',num2str(cfg.xlim(2))])
end 

posi        = get(subplot(2,5,10),'Position') % [left bottom width height]
c           = colorbar('Position', [posi(1)+pos(3)+0.135 posi(2)  0.03  posi(2)+posi(4)*2])
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'t-val')

end

end

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files']);
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);
rmpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);
