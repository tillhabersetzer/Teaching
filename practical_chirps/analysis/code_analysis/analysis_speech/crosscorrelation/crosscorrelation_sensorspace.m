close all; clear all; clc;

% Settings
%--------------------------------------------------------------------------
% choose subject 
% subjectlist = {'subject04'};
for i = 2:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjectlist{i-1} = [subject,num2str(i)]; 
end

% choose files
files2preproc = 'stories_maxfilter';

% choose envelope type for crosscorrelation
% envelopetype = 'onset_envelope';
envelopetype = {'onset_envelope','envelope'};
% envelopetype = {'envelope'};
% envelopetype = {'onset_envelope'};

% apply additional zscoring of all trials
zscoring = 1;

% load data with performed ica (1) or without (0)
ica_on = 1; 

% load data with whitening (1) or without (0)
% not necessary for correlation in sensorspace! cross correlation is
% normalized
whitening = 0; 

% bandpass neurodata
bp_freq = [0.5,45];

% downsampling frequency
fs_down = 250;

% lowpass audiodata
lp_freq_audio = 25;

% desired lags for crosscorrelation
timelags = [-100,900]; % in ms

% save data
save_data = 1;

% check data
check = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'))
% addpath for basic functions
addpath(fullfile('Z:','analysis','analysis_speech','helper_functions'));
% addpath for whitening functions
addpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));

%% preprocessing audiodata 
%--------------------------------------------------------------------------
audio_root        = fullfile('Z:','analysis','analysis_speech','audiodata');
path_audiodata{1} = fullfile(audio_root,'Das_Schwatzende_Herz','Das_schwatzende_Herz1_final.wav');
path_audiodata{2} = fullfile(audio_root,'Das_Schwatzende_Herz','Das_schwatzende_Herz2_final.wav');
path_audiodata{3} = fullfile(audio_root,'Die_Maske_des_Roten_Todes','Die_Maske_des_Roten_Todes1_final.wav');
path_audiodata{4} = fullfile(audio_root,'Die_Maske_des_Roten_Todes','Die_Maske_des_Roten_Todes2_final.wav');
 
% check if filter is okay
fsamp = 44100;
[b,a] = butter(3,lp_freq_audio/(fsamp/2),'low');
% h    = fvtool(b,a);
% h.Fs = fsamp;

for env = 1:length(envelopetype)
% hilbert/filter/derivative/half wave rectified
%----------------------------------------------
envelopes = cell(1,4);
for n = 1:4
    [raw_audiodata,fs] = audioread(path_audiodata{n}); 
    if ~isequal(fs,fsamp)
        error('unexpected sampling frequency!')
    end
    clear fs
    envelopes{n} = cal_envelope(raw_audiodata,fsamp,b,a,envelopetype{env});
    clear raw_audiodata
end

%% analysis of neurodata
%--------------------------------------------------------------------------
N_subj = length(subjectlist);

for s = 1:N_subj
eval(subjectlist{s})
% load data
%----------
if ica_on
    filename = [subjectdata.subjectname '_' files2preproc '_preproc_with_ica'];
else
    filename = [subjectdata.subjectname '_' files2preproc '_preproc_without_ica'];
end

if whitening
    filename = [filename '_whitened'];
end

data_preprocessed = importdata(fullfile(subjectdata.preprocessed_data_stories,[filename '.mat']));
info              = data_preprocessed.info;
data              = data_preprocessed.data_preprocessed;
clear data_preprocessed
F                 = length(data);

% append data of all stories 
%---------------------------
% dr                = data{1}.hdr;
cfg                = [];
cfg.keepsampleinfo = 'no';
data               = ft_appenddata(cfg,data{:});
% data.hdr           = hdr;

% filter data
%------------
if info.bp_freq(2) > bp_freq(2)
    cfg            = [];
    cfg.channel    = 'meg';         % 'all'
    cfg.detrend    = 'yes';
    cfg.demean     = 'yes';         % default: complete trial
    cfg.continuous = 'yes';
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = bp_freq;       % e.g. [0.1 200];
    data           = ft_resampledata(cfg,data);

end

% downsample data
%----------------
if ~isequal(data.fsample,fs_down)
    cfg            = [];
    cfg.resamplefs = fs_down;
    cfg.detrend    = 'no';
    data           = ft_resampledata(cfg,data);
end

%% covariance matrix for data check
%----------------------------------
if check
    shift  = 10;
    winlen = 10;
    [data_segments,~] = segment_data(data.trial,winlen,shift,fs_down);
    [time_segments,~] = segment_data(data.time,winlen,shift,fs_down);

    data2                = data;
    data2.time           = repmat(time_segments(1),1,length(time_segments));
    data2.trial          = data_segments;
    cfg                  = [];
    cfg.covariance       = 'yes';
    cfg.channel          = 'meg';
    cfg.covariancewindow = 'all';
    avg                  = ft_timelockanalysis(cfg,data2);
    clear data2 time_segments data_segments shift winlen
    
    % data covariance matrix
    %-----------------------
    selmag  = ft_chantype(avg.label, 'megmag');
    selgrad = ft_chantype(avg.label, 'megplanar');
    C = avg.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    figure
    subplot(1,2,1)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('MEG sensor covariance matrix')
    
    % singular values
    %----------------
    [~,sv,~] = svd(avg.cov);
    subplot(1,2,2)
    plot(log10(diag(sv)),'o');
    grid on
    title('Singular values of a MEG sensor covariance matrix')
    sgtitle(subjectdata.subjectname)
    clear avg
end

%% epoching and downsampling audiodata
% check which files are available in neurodata
idx      = contains({'story1part1_tsss_mc','story1part2_tsss_mc','story2part1_tsss_mc','story2part2_tsss_mc'},...
                    info.filenames);
envelope = envelopes(idx);

for f = 1:F % files               
    envelope{f} = epoch_audiodata(envelope{f},fsamp,info.triallength,1,fs_down);
end

% combine data
audiodata = [];
for f = 1:F
    audiodata = [audiodata,envelope{f}.trials];
end
clear envelope

%% match length of audio and neurodata
Nneuro = length(data.trial);
Naudio = length(audiodata);

if ~isequal(Nneuro,Naudio)
    error('error: number of epochs are different in neurodata and audiodata')
end

check_length = zeros(Naudio,2);
for n = 1:Naudio
    check_length(n,1) = size(data.trial{n},2);
    check_length(n,2) = length(audiodata{n});
end

% get index of minimum length for each trial
[M,~] = min(check_length,[],2);

% match length
for n = 1:Naudio
    data.trial{n} = data.trial{n}(:,1:M(n));
    audiodata{n}  = audiodata{n}(1:M(n));
    data.time{n}  = data.time{n}(1:M(n));
end
%% further segmentation of temporal aligned data 
winlen = 10;
shift  = 3; % 3 -> ~18 GB, 2 -> ~26 GB
[data.trial,~] = segment_data(data.trial,winlen,shift,fs_down);
[data.time,~]  = segment_data(data.time,winlen,shift,fs_down);
[audiodata,~]  = segment_data(audiodata,winlen,shift,fs_down);

%% apply zscoring
if zscoring
    N_segments = length(data.trial);
    for z = 1:N_segments
        data.trial{z} = zscore(data.trial{z},0,2);
        audiodata{z}  = zscore(audiodata{z},0,2);
    end
end
%% calculate crosscorrelation function
N_channel       = length(data.label);
N_segments      = length(data.trial);
timelags_sample = round(timelags/1000*fs_down); % in samples

% get correct index for lag vector
[~,lags]     = xcorr(data.trial{1}(1,:),audiodata{1},timelags_sample(2),'coeff');
idx_lags     = dsearchn(lags',timelags_sample');
timeshift    = lags(idx_lags(1):idx_lags(2))*(1/fs_down); 
rng('shuffle')
idx_shuffle = randperm(N_segments); % for random mapping of audiodata

% check if no elements are identical, so the difference should never be 0
while ~all(idx_shuffle-(1:N_segments)) % check for nonzero elements
    rng('shuffle')
    idx_shuffle = randperm(N_segments);
end

data.time         = repmat({timeshift},1,N_segments);
data_shuffle      = rmfield(data,{'trial','time'});
data_shuffle.time = repmat({timeshift},1,N_segments);

for e = 1:N_segments
    cfc  = zeros(N_channel,length(timeshift));
    cfcs = zeros(size(cfc));
    for c = 1:N_channel
        % xcorr(x(n+m),y(n),maxlag), x arrives latter for m>0 (so neuro and
        % postiv m)
        [r,~]         = xcorr(data.trial{e}(c,:),audiodata{e},timelags_sample(2),'coeff');
        [r_shuffle,~] = xcorr(data.trial{e}(c,:),audiodata{idx_shuffle(e)},timelags_sample(2),'coeff');
 
        cfc(c,:)      = r(idx_lags(1):idx_lags(2));
        cfcs(c,:)     = r_shuffle(idx_lags(1):idx_lags(2));
        disp(['channel: ',num2str(c),' | segment: ' ,num2str(e)])
    end
    data.trial{e}         = cfc;
    data_shuffle.trial{e} = cfcs;   
end

%% averages
cfg         = [];
cfg.channel = 'all';
avg         = ft_timelockanalysis(cfg,data);
avg_shuffle = ft_timelockanalysis(cfg,data_shuffle);

if check
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg,avg,avg_shuffle);
    sgtitle('average: all channels') 
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306mag_helmet.mat';
    ft_multiplotER(cfg,avg,avg_shuffle);
    sgtitle('average: magnetometer') 
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306planar_helmet.mat';
    ft_multiplotER(cfg,avg,avg_shuffle);
    sgtitle('average: gradiometer') 
end

%% statistic
% % Remark: crosscorrelation normalizes all data channels (mag,grad). So even
% % when the data is not prewhitened before, it could be processed with all
% % it's channels in one chunk ?!
% cfg                  = [];
% cfg.channel          = 'megmag';
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_indepsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan        = 2;
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 100;
% cfg_neighb           = [];
% cfg_neighb.channel   = 'megmag';
% cfg_neighb.method    = 'distance';
% cfg.neighbours       = ft_prepare_neighbours(cfg_neighb,avg.grad);
% %ft_neighbourplot(cfg,datapre_w)
% 
% design                 = zeros(1,2*N_segments);
% design(1,1:N_segments) = 1;
% design(1,N_segments)   = 2;
% 
% cfg.design = design;
% cfg.ivar   = 1; % row of design matrix that contains independent variable 
% 
% stat = ft_timelockstatistics(cfg,data,data_shuffle);

%%

if save_data
    new_dir = fullfile(subjectdata.speech,'crosscorrelation',envelopetype{env},'sensorlevel');
    if ~exist(new_dir, 'dir')
       mkdir(new_dir)
    end

    if ica_on
           filename_new = [subjectdata.subjectname '_' files2preproc '_crosscorr_ica'];
       else
           filename_new = [subjectdata.subjectname '_' files2preproc '_crosscorr_noica'];
    end
    
    if whitening
        filename_new = [filename_new '_whitened'];
    end
    
    if zscoring
        filename_new = [filename_new '_zscored'];
    end
    save(fullfile(new_dir,[filename_new,'.mat']),'avg','avg_shuffle','info','-v7.3')  
end

end

end % envelopetype

%% Clean up
rmpath(fullfile('Z:','analysis','subject_files'))
rmpath(fullfile('Z:','analysis','analysis_speech','helper_functions'));
rmpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));
