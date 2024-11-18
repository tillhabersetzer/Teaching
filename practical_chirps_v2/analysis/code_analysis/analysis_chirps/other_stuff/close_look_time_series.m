close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For statistical testing between activation vs baseline 
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
bp_freq = [0.1,120];

% apply ica
ica_on = 1;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% downsample data
downsample_data = 0;

% show plot
show_plots = 0;

% choose statistic method 
cmb_plr = 1; % combined planar gradiometers (noise sensitive!!!)

plr_1   = 1; % normal planar gradiometers 

plr_2   = 1; % planar gradiometers, but baseline trials are replaced by averaged baseline
             % so activation is checked againt time averaged basline
             
mag     = 1; % magnetometers, but baseline trials are replaced by averaged baseline
             % so activation is checked againt time averaged basline

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
    cfg.trialdef.prestim    = 0.55;                  % in seconds
    cfg.trialdef.poststim   = 0.55;                  % in seconds
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
end

new_dir = [subjectdata.chirps_statistics filesep files2preproc];
if ~exist(new_dir, 'dir')
    mkdir(new_dir)
end

%% planar gradiometers statistics I
if plr_1
%  compute the statistical evaluation using permutation approach - within 
cfg        = [];
cfg.toilim = [-0.5 0];
datapre    = ft_redefinetrial(cfg,data_preprocessed);
cfg.toilim = [0 0.5];
datapst    = ft_redefinetrial(cfg,data_preprocessed);

datapre.time = datapst.time;

% statistic
cfg                  = [];
cfg.channel          = 'megplanar';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb           = [];
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'neuromag306planar_neighb.mat';
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

savefig([new_dir filesep 'raweffect_megplanar_1.fig'])

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
   % pos_int(i1) = any(pos(i2, m(k):m(k+1)), 2);
   % neg_int(i1) = any(neg(i2, m(k):m(k+1)), 2);

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
savefig([new_dir filesep 'cluster_permutation_megplanar_1.fig'])

end

%% planar gradiometers statistics II
if plr_2
% compute the statistical evaluation using permutation approach - within 
%% my idea......
% The statistic compares the power in every sample (i.e., a (channel,time)-triplet) 
% in the activation period with the corresponding time-averaged power 
% (i.e., the average over the temporal dimension) in the baseline period. 
% i try this because cfg.statistic =  'ft_statfun_actvsblT' seems to work
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
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb           = [];
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'neuromag306planar_neighb.mat';
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

savefig([new_dir filesep 'raweffect_megplanar_2.fig'])

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

savefig([new_dir filesep 'cluster_permutation_megplanar_2.fig'])

end

%% combined planar gradiometers statistics 
if cmb_plr

%!!! maybe not the best way to work with single combined gradiometer
% channels because they contain a lot of noise and they are all positive, so
% they don't average out!!!

% compute the statistical evaluation using permutation approach - within 

% demean data_preprocessed to reduce summation of noise (x^2+y^2)
cfg                      = [];
cfg.demean               = 'yes';
data_preprocessed_demean = ft_preprocessing(cfg,data_preprocessed);

cfg                   = [];
cfg.method            = 'sum';
cfg.updatesens        = 'yes';        
data_preprocessed_cmb = ft_combineplanar(cfg,data_preprocessed_demean);

%demean signals (pre+pst together)
cfg                   = [];
cfg.demean            = 'yes';
data_preprocessed_cmb = ft_preprocessing(cfg,data_preprocessed_cmb);

% seperate data
cfg        = [];
cfg.toilim = [-0.5 0];
datapre    = ft_redefinetrial(cfg,data_preprocessed_cmb);
cfg.toilim = [0 0.5];
datapst    = ft_redefinetrial(cfg,data_preprocessed_cmb);

datapre.time = datapst.time;

% optional - check against time averaged baseline
datapre.trial = cellfun(@(x) repmat(mean(x,2),1,size(x,2)),datapre.trial,'UniformOutput',false);

% statistic
cfg                  = [];
cfg.channel          = 'megplanar';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb           = [];
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'neuromag306cmb_neighb.mat';
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

savefig([new_dir filesep 'raweffect_combined_megplanar.fig'])

% have a look at raweffect
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg, raweffect);
suptitle(['raweffect cmb'])

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
   % pos_int     = zeros(numel(raweffect.label),1);
   % neg_int     = zeros(numel(raweffect.label),1);
   % pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   % neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
   % less severe
   pos_int(i1) = any(pos(i2, m(k):m(k+1)), 2);
   neg_int(i1) = any(neg(i2, m(k):m(k+1)), 2);

   cfg.highlight        = 'on';
   % Get the index of each significant channel
   cfg.highlightchannel = find(pos_int | neg_int);
   %cfg.highlightchannel = find(pos_int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 10;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306cmb.lay';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg, raweffect);
end

savefig([new_dir filesep 'cluster_permutation_combined_megplanar.fig'])

% have a look at datapre and datapst
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg, datapre);
suptitle([''])

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg, datapst);
suptitle([''])

end

%% magnetometers statistics II
if mag
% compute the statistical evaluation using permutation approach - within 
%% my idea......
% The statistic compares the power in every sample (i.e., a (channel,time)-triplet) 
% in the activation period with the corresponding time-averaged power 
% (i.e., the average over the temporal dimension) in the baseline period. 
% i try this because cfg.statistic = 'ft_statfun_actvsblT' seems to work
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
cfg.channel          = 'megmag';
%cfg.latency          = [0 0.5];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; % dependent samples T-statistic 
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb           = [];
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'neuromag306mag_neighb.mat';
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
cfg.channel = 'megmag';
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
cfg.layout     = 'neuromag306mag.lay';
ft_multiplotER(cfg, raweffect);
suptitle(['raweffect'])

savefig([new_dir filesep 'raweffect_megmag.fig'])

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
   cfg.layout           = 'neuromag306mag.lay';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg, raweffect);
end

savefig([new_dir filesep 'cluster_permutation_megmag.fig'])

end

close all

end

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])

