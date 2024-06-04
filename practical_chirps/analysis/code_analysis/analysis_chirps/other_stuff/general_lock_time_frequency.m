close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the avarages and
% beamforming
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
%bp_freq = [4, 120];

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
    cfg.trialdef.prestim    = 0.5;                  % in seconds
    cfg.trialdef.poststim   = 0.5;                   % in seconds
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
    %cfg.baselinewindow = [-0.5 0];
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
    %cfg.bpfilter  = 'yes';
    %cfg.bpfreq    = bp_freq;      
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
cfg            = [];
cfg.method     = 'sum';
cfg.updatesens = 'yes'; % need old information for dipole fitting 'no'
data_avg_cmb   = ft_combineplanar(cfg,data_avg);

%% just to test - combine gradiometers before averaging
cfg                   = [];
cfg.method            = 'sum';
cfg.updatesens        = 'yes'; 
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
    cfg.xlim       = [-0.5, 0.5];
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

%% Computing time-frequency power representations

% wavelet
%--------
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'wavelet';
cfg.width      = 5;
cfg.foi        = 3:2:30;
cfg.toi        = -0.5:.05:0.5;
cfg.keeptrials = 'no';
tfr_wavelet    = ft_freqanalysis(cfg,  data_preprocessed);
% baseline correction
cfg              =[];
cfg.baseline     =[-0.5 -0.1];
cfg.baselinetype = 'absolute';
tfr_wavelet_bl   = ft_freqbaseline(cfg, tfr_wavelet);

% multitaper time-frequency transformation
%-----------------------------------------
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = 3:3:30;
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.34;
cfg.toi        = -0.5:.05:0.5;
cfg.keeptrials = 'no';
tfr            = ft_freqanalysis(cfg,  data_preprocessed);
% baseline correction
cfg              =[];
cfg.baseline     =[-0.5 -0.1];
cfg.baselinetype = 'absolute';
tfr_bl           = ft_freqbaseline(cfg, tfr);

%% for averaged combined gradiometers
% wavelet
%--------
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'wavelet';
cfg.width      = 5;
cfg.foi        = 3:2:30;
cfg.toi        = -0.5:.05:0.5;
cfg.keeptrials = 'no';
tfr_wavelet_avg_cmb = ft_freqanalysis(cfg,  data_avg_cmb);
% baseline correction
cfg              =[];
cfg.baseline     =[-0.5 -0.1];
cfg.baselinetype = 'absolute';
tfr_wavelet_bl_avg_cmb   = ft_freqbaseline(cfg, tfr_wavelet_avg_cmb);

% multitaper time-frequency transformation
%-----------------------------------------
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = 3:3:30;
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.34;
cfg.toi        = -0.5:.05:0.5;
cfg.keeptrials = 'no';
tfr_avg_cmb    = ft_freqanalysis(cfg,  data_avg_cmb);
% baseline correction
cfg              =[];
cfg.baseline     =[-0.5 -0.1];
cfg.baselinetype = 'absolute';
tfr_bl_avg_cmb   = ft_freqbaseline(cfg, tfr_avg_cmb);

%% Plot
if show_plots  
    
    choice = 'tfr'; % 'tfr'
    if strcmp(choice,'wavelet')
        tfr2plot = tfr_wavelet_bl
        tfr2plot_avg_cmb = tfr_wavelet_bl_avg_cmb
    elseif strcmp(choice,'tfr')
        tfr2plot = tfr_bl
        tfr2plot_avg_cmb = tfr_bl_avg_cmb
    end
    % magnetometer
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306mag.lay';
    cfg.xlim       = [-0.5, 0.5];
    ft_multiplotTFR(cfg, tfr2plot)
    suptitle([choice ': average magnetometer'])

    figure
    cfg          = [];
    cfg.xlim     = [0.07, 0.140];
    cfg.ylim     = [4 15]; % frequency lim
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306mag.lay';
    ft_topoplotTFR(cfg, tfr2plot);
    suptitle([choice ': topoplot magnetometer'])

    figure
    cfg        = [];
    cfg.xlim   = [-0.1 : 0.1 : 0.3];  % Define time intervals
    cfg.ylim   = [4 15]; % frequency lim
    cfg.layout = 'neuromag306mag.lay';
    ft_topoplotTFR(cfg, tfr2plot);
    suptitle([choice ': topoplot series magnetometer'])
    
    % gradiometer
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306planar.lay';
    cfg.xlim       = [-0.5, 0.5];
    ft_multiplotTFR(cfg, tfr2plot)
    suptitle([choice ': average gradiometer'])

    figure
    cfg          = [];
    cfg.xlim     = [0.07, 0.140];
    cfg.ylim     = [4 15]; % frequency lim
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306planar.lay';
    ft_topoplotTFR(cfg, tfr2plot);
    suptitle([choice ': topoplot gradiometer'])

    figure
    cfg        = [];
    cfg.xlim   = [-0.1 : 0.1 : 0.3];  % Define time intervals
    cfg.ylim   = [4 15]; % frequency lim
    cfg.layout = 'neuromag306planar.lay';
    ft_topoplotTFR(cfg, tfr2plot);
    suptitle([choice ': topoplot series gradiometer'])
    
    % averaged combined gradiometer
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb.lay';
    cfg.xlim       = [-0.5, 0.5];
    ft_multiplotTFR(cfg, tfr2plot_avg_cmb)
    suptitle([choice ': averaged combined gradiometer'])

    figure
    cfg          = [];
    cfg.xlim     = [0.07, 0.140];
    cfg.ylim     = [4 15]; % frequency lim
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306cmb.lay';
    ft_topoplotTFR(cfg, tfr2plot_avg_cmb);
    suptitle([choice ': topoplot averaged combined gradiometer'])

    figure
    cfg        = [];
    cfg.xlim   = [-0.1 : 0.1 : 0.3];  % Define time intervals
    cfg.ylim   = [4 15]; % frequency lim
    cfg.layout = 'neuromag306cmb.lay';
    ft_topoplotTFR(cfg, tfr2plot_avg_cmb);
    suptitle([choice ': topoplot series averaged combined gradiometer'])
    
    % plot also a single channel
    figure
    cfg = [];
    cfg.channel = 'MEG1431';
    subplot(2,2,1)
    ft_singleplotTFR(cfg, tfr2plot);
    cfg.channel = 'MEG1332+1333';
    subplot(2,2,2)
    ft_singleplotTFR(cfg, tfr2plot_avg_cmb);
    cfg.channel = 'MEG2422';
    subplot(2,2,3)
    ft_singleplotTFR(cfg, tfr2plot);
    cfg.channel = 'MEG2423';
    subplot(2,2,4)
    ft_singleplotTFR(cfg, tfr2plot);
    suptitle([choice ': selected channels'])
    
end

%% try another approach
% combine gradiometers after freqanalysis
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = 3:3:30;
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.34;
cfg.toi        = -0.5:.05:0.5;
cfg.keeptrials = 'yes';
tfr1           = ft_freqanalysis(cfg, data_preprocessed);
% baseline correction
cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
tfr1_bl          = ft_freqbaseline(cfg, tfr1);

% combine gradiometers
cfg            = [];
cfg.method     = 'sum';
cfg.updatesens = 'no';
tfr_before_cmb = ft_combineplanar(cfg,tfr1_bl);

% average for raweffect
cfg            = [];
tfr_before_cmb = ft_freqdescriptives(cfg, tfr_before_cmb);

if show_plots
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb.lay';
    cfg.xlim       = [-0.5, 0.5];
    ft_multiplotTFR(cfg, tfr_before_cmb)
    suptitle(['tfr->cmb->avg'])

    figure
    cfg          = [];
    cfg.xlim     = [0.07, 0.140];
    cfg.ylim     = [4 15]; % frequency lim
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306cmb.lay';
    ft_topoplotTFR(cfg, tfr_before_cmb);
    suptitle(['tfr->cmb->avg'])

    figure
    cfg        = [];
    cfg.xlim   = [-0.1 : 0.1 : 0.3];  % Define time intervals
    cfg.ylim   = [4 15]; % frequency lim
    cfg.layout = 'neuromag306cmb.lay';
    ft_topoplotTFR(cfg, tfr_before_cmb);
    suptitle(['tfr->cmb->avg'])  
end



%% Within subject statistics on time-frequency representations / Within trial statictics
% Sometimes one might be interested in power modulations significantly 
% different from pre stimulus baseline. In the following section we will illustrate how to achieve this.

cfg        = [];
cfg.toilim = [-0.5 0];
datapre    = ft_redefinetrial(cfg,data_preprocessed);
cfg.toilim = [0 0.5];
datapst    = ft_redefinetrial(cfg,data_preprocessed);
%%
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = 4:2:30;
%cfg.method     = 'wavelet';
%cfg.width      = 3;
cfg.toi        = -0.5:.05:0.5;
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.25;
cfg.toi        = -0.5:.05:0;
cfg.keeptrials = 'yes';
tfrpre         = ft_freqanalysis(cfg, datapre);
cfg.toi        = 0:.05:0.5;
tfrpst         = ft_freqanalysis(cfg, datapst);

% combine gradiometers
cfg            = [];
cfg.method     = 'sum';
cfg.updatesens = 'no';
tfrpre         = ft_combineplanar(cfg,tfrpre);
tfrpst         = ft_combineplanar(cfg,tfrpst);

tfrpre.time = tfrpst.time;
tfrpre.freq = round(tfrpre.freq);
tfrpst.freq = round(tfrpst.freq);

if show_plots
    % average for raweffect
    cfg        = [];
    tfrpst_avg = ft_freqdescriptives(cfg, tfrpst);

    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb.lay';
    cfg.xlim       = [-0.5, 0.5];
    ft_multiplotTFR(cfg, tfrpst_avg)
    suptitle([''])

    figure
    cfg          = [];
    cfg.xlim     = [0.07, 0.140];
    cfg.ylim     = [4 15]; % frequency lim
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306cmb.lay';
    ft_topoplotTFR(cfg, tfrpst_avg);
    suptitle([''])

    figure
    cfg        = [];
    %cfg.xlim   = [-0.1 : 0.1 : 0.3];  % Define time intervals
    cfg.ylim   = [4 15]; % frequency lim
    cfg.layout = 'neuromag306cmb.lay';
    ft_topoplotTFR(cfg, tfrpst_avg);
    suptitle(['']) 
end

% Now we compute the statistical evaluation using permutation approach
cfg = [];
cfg.channel          = 'meggrad';
%cfg.latency          = [0 0.5];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_actvsblT'; % activation-versus-baseline T-statistic
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
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, tfrpre);

ntrials                       = size(tfrpst.powspctrm,1);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable
cfg.uvar   = 2; % row of design matrix that contains unit of observation

stat = ft_freqstatistics(cfg, tfrpst, tfrpre);

cfg               = [];
cfg.channel       = {'MEG2522+2523'};
%cfg.channel       = 'MEG2423';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle     = 'outline';
%cfg.zlim          = [-5 5];
%cfg.xlim          = [0 0.5];
figure;
subplot(1,2,1)
ft_singleplotTFR(cfg,stat);
cfg           = [];
cfg.layout    = 'neuromag306cmb.lay';
cfg.parameter = 'stat';
%cfg.xlim      = [0.06 0.3];
%cfg.ylim      = [5 8];
%cfg.zlim      = [-5 5]
subplot(1,2,2)
ft_topoplotTFR(cfg,stat);














end

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])




