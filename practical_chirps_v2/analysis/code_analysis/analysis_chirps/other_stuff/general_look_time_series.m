close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For statistical testing between activation vs baseline 
% NO baseline correction should be applied to the data
%
% Inspired by: http://www.fieldtriptoolbox.org/tutorial/salzburg/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject04'};

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [4, 120];

% apply ica
ica_on = 1;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% downsample data
downsample_data = 0;

% show plot
show_plots = 0;
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

    %% define trials
    cfg                     = [];
    cfg.dataset             = path_dataset;
    cfg.trialfun            = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = 2;
    cfg.trialdef.prestim    = 0.55;                  % in seconds
    cfg.trialdef.poststim   = 0.55;                   % in seconds
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
    
    %% preprocess data 
    cfg                = [];
    cfg.dataset        = path_dataset;
    cfg.trl            = trl_new;
    cfg.channel        = 'meg'; 
    cfg.coilaccuracy   = 0;           % ensure that sensors are expressed in SI units
    %cfg.baselinewindow = [-0.55 0];
    %cfg.demean         = 'yes';          
    data               = ft_preprocessing(cfg);
    
    %% reject earlier specified independet components
    if ica_on
        data = reject_independent_components(data,subjectdata,ica_files);
    end
       
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
    
    %% filter data
    cfg           = [];
    %cfg.padtype   = 'mirror';
    %cfg.padding    = 1;
    cfg.channel   = 'meg';        
    cfg.bpfilter  = 'yes';
    cfg.bpfreq    = bp_freq;      
    cfg.dftfilter = 'yes';        
    cfg.dftfreq   = [50 100 150];  
    data          = ft_preprocessing(cfg, data);
    
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

%% average data
cfg      = [];
data_avg = ft_timelockanalysis(cfg,data_preprocessed);

%% combine planar gradient
cfg                = [];
cfg.method         = 'sum';
cfg.updatesens     = 'yes'; % need old information for dipole fitting 'no'
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.5 0];
data_avg_cmb       = ft_combineplanar(cfg,data_avg);

%% just to test - combine gradiometers before averaging
cfg                   = [];
cfg.method            = 'sum';
cfg.updatesens        = 'yes'; 
cfg.demean            = 'yes';
cfg.baselinewindow    = [-0.4 0];
data_preprocessed_cmb = ft_combineplanar(cfg,data_preprocessed);

cfg          = [];
data_cmb_avg = ft_timelockanalysis(cfg,data_preprocessed_cmb);

%% Plot
if show_plots
    % magnetometer
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306mag.lay';
    %cfg.xlim       = [-0.55, 0.55];
    ft_multiplotER(cfg,data_avg_cmb);
    suptitle('avergage magnetometer')

    figure
    cfg          = [];
    cfg.xlim     = [0.07, 0.140];
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306mag.lay';
    ft_topoplotER(cfg,data_avg_cmb);
    suptitle('topoplot magnetometer')

    figure
    cfg        = [];
    cfg.xlim   = [-0.1 : 0.1 : 0.5];  % Define time intervals
    cfg.layout = 'neuromag306mag.lay';
    ft_topoplotER(cfg,data_avg_cmb);
    suptitle('topoplot series magnetometer')

    % planar gradient - averaged - after combining
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb.lay';
    ft_multiplotER(cfg, data_avg_cmb);
    suptitle('average combined gradiometers after avg')

    figure
    cfg          = [];
    cfg.xlim     = [0.07, 0.140];
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306cmb.lay';
    ft_topoplotER(cfg,data_avg_cmb);
    suptitle('topoplot combined gradiometers after avg')

    figure
    cfg        = [];
    cfg.xlim   = [-0.1 : 0.1 : 0.5];  % Define time intervals
    cfg.layout = 'neuromag306cmb.lay';
    ft_topoplotER(cfg,data_avg_cmb);
    suptitle('topoplot series combined gradiometers after avg')
    
    % planar gradient - combined before averaging
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb.lay';
    ft_multiplotER(cfg, data_cmb_avg);
    suptitle('average combined gradiometers before avg')

    figure
    cfg          = [];
    cfg.xlim     = [0.07, 0.140];
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306cmb.lay';
    ft_topoplotER(cfg,data_cmb_avg);
    suptitle('topoplot combined gradiometers before avg')

    figure
    cfg        = [];
    cfg.xlim   = [-0.1 : 0.1 : 0.5];  % Define time intervals
    cfg.layout = 'neuromag306cmb.lay';
    ft_topoplotER(cfg,data_cmb_avg);
    suptitle('topoplot series combined gradiometers before avg')
end

%% planar gradiometers
%% compute the statistical evaluation using permutation approach - between
cfg        = [];
cfg.toilim = [-0.5 0];
datapre    = ft_redefinetrial(cfg,data_preprocessed);
cfg.toilim = [0 0.5];
datapst    = ft_redefinetrial(cfg,data_preprocessed);

datapre.time = datapst.time;

% statistic
cfg                  = [];
cfg.channel          = 'megplanar';
%cfg.latency          = [0 0.5];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT'; % independent samples T-statistic 
% This statistic compares the power in every sample (i.e., a (channel,frequency,time)-triplet) 
% in the activation period with the corresponding time-averaged power (i.e., the average over the temporal dimension)
% in the baseline period. The comparison of the activation and the time-averaged baseline power is performed 
% by means of a dependent samples T-statistic.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb           = [];
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'neuromag306planar_neighb.mat';
%cfg_neighb.template  = 'neuromag306cmb_neighb.mat';
%cfg_neighb.template  = 'neuromag306mag_neighb.mat';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, datapre);

ntrials                       = length(datapre.trial);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable

stat = ft_timelockstatistics(cfg, datapst, datapre);

% raw effect
cfg         = [];
cfg.channel = 'megplanar';
datapre_avg = ft_timelockanalysis(cfg, datapre);
datapst_avg = ft_timelockanalysis(cfg, datapst);

cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffect     = ft_math(cfg,datapst_avg, datapre_avg);

% have a look at raweffect
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
ft_multiplotER(cfg, raweffect);
suptitle(['raweffect'])

% make nice plots
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha); % use all significant clusters
pos              = ismember(stat.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < stat.cfg.alpha);
neg               = ismember(stat.negclusterslabelmat, neg_signif_clust);

timestep      = 0.025; 
sampling_rate = datapre.fsample; 
sample_count  = length(stat.time);
j = [0:timestep:0.5]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in samples

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect.label, stat.label);

figure
for k = 1:20
   subplot(4,5,k);
   cfg      = [];
   cfg.xlim = [j(k) j(k+1)];  
   % cfg.zlim = [-2.5e-13 2.5e-13];
   % Next, check which channels are significant over the
   % entire time interval of interest.
   pos_int     = zeros(numel(raweffect.label),1);
   neg_int     = zeros(numel(raweffect.label),1);
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
   cfg.layout           = 'neuromag306planar.lay';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg, raweffect);
end

%% planar gradiometers
%% compute the statistical evaluation using permutation approach - within 

cfg        = [];
cfg.toilim = [-0.5 0];
datapre    = ft_redefinetrial(cfg,data_preprocessed);
cfg.toilim = [0 0.5];
datapst    = ft_redefinetrial(cfg,data_preprocessed);

datapre.time = datapst.time;

% statistic
cfg                  = [];
cfg.channel          = 'megplanar';
%cfg.latency          = [0 0.5];
cfg.method           = 'montecarlo';
%cfg.statistic        = 'ft_statfun_actvsblT';
cfg.statistic        = 'ft_statfun_depsamplesT'; % dependent samples T-statistic 
% This statistic compares the power in every sample (i.e., a (channel,frequency,time)-triplet) 
% in the activation period with the corresponding time-averaged power (i.e., the average over the temporal dimension)
% in the baseline period. The comparison of the activation and the time-averaged baseline power is performed 
% by means of a dependent samples T-statistic.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb           = [];
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'neuromag306planar_neighb.mat';
%cfg_neighb.template  = 'neuromag306cmb_neighb.mat';
%cfg_neighb.template  = 'neuromag306mag_neighb.mat';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, datapre);

ntrials                       = length(datapre.trial);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable
cfg.uvar   = 2; % row of design matrix that contains unit of observation

stat = ft_timelockstatistics(cfg, datapst, datapre);

% raw effect
cfg         = [];
cfg.channel = 'megplanar';
datapre_avg = ft_timelockanalysis(cfg, datapre);
datapst_avg = ft_timelockanalysis(cfg, datapst);

cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffect     = ft_math(cfg,datapst_avg, datapre_avg);

% have a look at raweffect
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
ft_multiplotER(cfg, raweffect);
suptitle(['raweffect'])

% make nice plots
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha); % use all significant clusters
pos              = ismember(stat.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < stat.cfg.alpha);
neg               = ismember(stat.negclusterslabelmat, neg_signif_clust);

timestep      = 0.025; 
sampling_rate = datapre.fsample; 
sample_count  = length(stat.time);
j = [0:timestep:0.5]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in samples

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect.label, stat.label);

figure
for k = 1:20
   subplot(4,5,k);
   cfg      = [];
   cfg.xlim = [j(k) j(k+1)];  
   % cfg.zlim = [-2.5e-13 2.5e-13];
   % Next, check which channels are significant over the
   % entire time interval of interest.
   pos_int     = zeros(numel(raweffect.label),1);
   neg_int     = zeros(numel(raweffect.label),1);
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
   cfg.layout           = 'neuromag306planar.lay';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg, raweffect);
end

%% planar gradiometers
%% compute the statistical evaluation using permutation approach - within 
%% my idea......
% The statistic compares the power in every sample (i.e., a (channel,time)-triplet) 
% in the activation period with the corresponding time-averaged power 
% (i.e., the average over the temporal dimension) in the baseline period. 
% i try this because cfg.statistic = 'ft_statfun_depsamplesT'seems to work
% only for time-frequency data in ft_freqstatistics

cfg        = [];
cfg.toilim = [-0.5 0];
datapre    = ft_redefinetrial(cfg,data_preprocessed);
cfg.toilim = [0 0.5];
datapst    = ft_redefinetrial(cfg,data_preprocessed);

% so replace datapre with time averaged data
% average over time (columns) and replace matrix data with average
datapre.trial = cellfun(@(x) repmat(mean(x,2),1,size(x,2)),datapre.trial,'UniformOutput',false);

datapre.time = datapst.time;

% statistic
cfg                  = [];
cfg.channel          = 'megplanar';
%cfg.latency          = [0 0.5];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; % dependent samples T-statistic 
% This statistic compares the power in every sample (i.e., a (channel,frequency,time)-triplet) 
% in the activation period with the corresponding time-averaged power (i.e., the average over the temporal dimension)
% in the baseline period. The comparison of the activation and the time-averaged baseline power is performed 
% by means of a dependent samples T-statistic.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb           = [];
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'neuromag306planar_neighb.mat';
%cfg_neighb.template  = 'neuromag306cmb_neighb.mat';
%cfg_neighb.template  = 'neuromag306mag_neighb.mat';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, datapre);

ntrials                       = length(datapre.trial);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable
cfg.uvar   = 2; % row of design matrix that contains unit of observation

stat = ft_timelockstatistics(cfg, datapst, datapre);

% raw effect
cfg         = [];
cfg.channel = 'megplanar';
datapre_avg = ft_timelockanalysis(cfg, datapre);
datapst_avg = ft_timelockanalysis(cfg, datapst);

cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffect     = ft_math(cfg,datapst_avg, datapre_avg);

% have a look at raweffect
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
ft_multiplotER(cfg, raweffect);
suptitle(['raweffect'])

% make nice plots
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha); % use all significant clusters
pos              = ismember(stat.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < stat.cfg.alpha);
neg               = ismember(stat.negclusterslabelmat, neg_signif_clust);

timestep      = 0.025; 
sampling_rate = datapre.fsample; 
sample_count  = length(stat.time);
j = [0:timestep:0.5]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in samples

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect.label, stat.label);

figure
for k = 1:20
   subplot(4,5,k);
   cfg      = [];
   cfg.xlim = [j(k) j(k+1)];  
   % cfg.zlim = [-2.5e-13 2.5e-13];
   % Next, check which channels are significant over the
   % entire time interval of interest.
   pos_int     = zeros(numel(raweffect.label),1);
   neg_int     = zeros(numel(raweffect.label),1);
   pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
   % less severe
   %pos_int(i1) = any(pos(i2, m(k):m(k+1)), 2);
   %neg_int(i1) = any(neg(i2, m(k):m(k+1)), 2);

   cfg.highlight        = 'on';
   % Get the index of each significant channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 10;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306planar.lay';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg, raweffect);
end

%% combined planar gradiometers
%% compute the statistical evaluation using permutation approach - within 

cfg        = [];
cfg.toilim = [-0.5 0];
datapre    = ft_redefinetrial(cfg,data_preprocessed);
cfg.toilim = [0 0.5];
datapst    = ft_redefinetrial(cfg,data_preprocessed);

% combine planar gradient
cfg            = [];
cfg.method     = 'sum';
cfg.updatesens = 'yes'; % need old information for dipole fitting 'no'
datapre_cmb    = ft_combineplanar(cfg,datapre);
datapst_cmb    = ft_combineplanar(cfg,datapst);

datapre_cmb.time = datapst_cmb.time;

% statistic
cfg                  = [];
cfg.channel          = 'megplanar';
%cfg.latency          = [0 0.5];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; % dependent samples T-statistic 
% This statistic compares the power in every sample (i.e., a (channel,frequency,time)-triplet) 
% in the activation period with the corresponding time-averaged power (i.e., the average over the temporal dimension)
% in the baseline period. The comparison of the activation and the time-averaged baseline power is performed 
% by means of a dependent samples T-statistic.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb           = [];
cfg_neighb.method    = 'template';
%cfg_neighb.template  = 'neuromag306planar_neighb.mat';
cfg_neighb.template  = 'neuromag306cmb_neighb.mat';
%cfg_neighb.template  = 'neuromag306mag_neighb.mat';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, datapre_cmb);

ntrials                       = length(datapre_cmb.trial);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable
cfg.uvar   = 2; % row of design matrix that contains unit of observation

stat = ft_timelockstatistics(cfg, datapst_cmb, datapre_cmb);

% raw effect
cfg             = [];
cfg.channel     = 'megplanar';
datapre_cmb_avg = ft_timelockanalysis(cfg, datapre_cmb);
datapst_cmb_avg = ft_timelockanalysis(cfg, datapst_cmb);

cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffect_cmb = ft_math(cfg,datapst_cmb_avg, datapre_cmb_avg);

% have a look at raweffect
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg, raweffect_cmb);
suptitle(['raweffect_cmb'])

% make nice plots
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha); % use all significant clusters
pos              = ismember(stat.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < stat.cfg.alpha);
neg               = ismember(stat.negclusterslabelmat, neg_signif_clust);

timestep      = 0.025; 
sampling_rate = datapre.fsample; 
sample_count  = length(stat.time);
j = [0:timestep:0.5]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in samples

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect_cmb.label, stat.label);

figure
for k = 1:20
   subplot(4,5,k);
   cfg      = [];
   cfg.xlim = [j(k) j(k+1)];  
   % cfg.zlim = [-2.5e-13 2.5e-13];
   % Next, check which channels are significant over the
   % entire time interval of interest.
   pos_int     = zeros(numel(raweffect_cmb.label),1);
   neg_int     = zeros(numel(raweffect_cmb.label),1);
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
   cfg.layout           = 'neuromag306cmb.lay';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg, raweffect_cmb);
end





end

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])




