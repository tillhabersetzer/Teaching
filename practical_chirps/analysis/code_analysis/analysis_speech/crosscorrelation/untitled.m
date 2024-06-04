close all; clear all; clc;

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjectlist = {'subject04'};

% choose files
files2preproc = 'stories_maxfilter';

% load data with performed ica (1) or without (0)
ica_status = 0; 

% load data with whitening (1) or without (0)
whitening = 0;

% downsampling frequency
fs_down = 250;

% lowpass audiodata
lp_freq_audio = 25;

% bandpass neurodata
bp_freq = [0.5,45];

% filewise zcore 
filewise_zscore = 0;

% trialwise_zcore 
trialwise_zscore = 0;

% shift length and window length for segmentation of data
shift  = 1;  % in sec
winlen = 10; % in sec

% save data
save_data = 1;

% check data
check = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'))
% addpath for basic functions
addpath(fullfile('Z:','analysis','preprocessing_batch','helper_functions'))
addpath(fullfile('Z:','analysis','analysis_speech','helper_functions'));

if filewise_zscore && trialwise_zscore
    error(['decide for one methode of z-scoring'])
end

%% preprocessing audiodata 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
audio_root        = fullfile('Z:','analysis','analysis_speech','audiodata');
path_audiodata{1} = fullfile(audio_root,'Das_Schwatzende_Herz','Das_schwatzende_Herz1_final.wav');
path_audiodata{2} = fullfile(audio_root,'Das_Schwatzende_Herz','Das_schwatzende_Herz2_final.wav');
path_audiodata{3} = fullfile(audio_root,'Die_Maske_des_Roten_Todes','Die_Maske_des_Roten_Todes1_final.wav');
path_audiodata{4} = fullfile(audio_root,'Die_Maske_des_Roten_Todes','Die_Maske_des_Roten_Todes2_final.wav');

% number of stories
F = 4;
 
% check if filter is okay
fsamp = 44100;
[b,a] = butter(3,lp_freq_audio/(fsamp/2),'low');
% h    = fvtool(b,a);
% h.Fs = fsamp;

onset_envelope  = cell(1,F);
raw_audiodata   = cell(1,F);

%% hilbert/filter/derivative/half wave rectified
for n = 1:F
    [raw_audiodata{n},fs] = audioread(path_audiodata{n}); 
    if ~isequal(fs,fsamp)
        error('unexpected sampling frequency!')
    end
    clear fs
    onset_envelope{n}    = cal_envelope(raw_audiodata{n},fsamp,b,a);
end

% have a look at speech onset envelope
%--------------------------------------------------------------------------
if check
    label = {'Das_schwatzende_Herz1','Das_schwatzende_Herz2','Die_Maske_des_Roten_Todes1','Die_Maske_des_Roten_Todes2'};
    % choose file (n=1,2,3)
    n = 1;
    x = (1:length(raw_audiodata{n}))*(1/fsamp);
    % dx is the width of the axis 'widxow'
    dx = 20; % in sec
    % Initialise the figure once, and we only need to set the properties once
    fig{n} = figure(1); clf;
    set( fig{n}, 'doublebuffer', 'on'); 
    % Create a placeholder for axes objects
    ax = gobjects(2,1);
    % Create plots, storing them in the axes object
    ax(1) = subplot(2,1,1);
    plot(x, raw_audiodata{n});
    title('raw audiodata')
    ax(2) = subplot(2,1,2);
    plot(x, onset_envelope{n}); 
    xlabel('t / s')
    title('speech onset envelope')
    % Set up the scroller for the array of axes objects in 'ax'
    scrollplot( dx, x, ax)
    sgtitle(label{n}) 
end
%--------------------------------------------------------------------------

%% preprocess neurodata 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over subjects 
%-------------------
N_subj  = length(subjectlist);

for s = 1:N_subj
eval(subjectlist{s})
% load data
%----------
if ica_status
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

% append data of all stories 
%---------------------------
% hdr                = data{1}.hdr;
cfg                = [];
cfg.keepsampleinfo = 'no';
data               = ft_appenddata(cfg,data{:});
% data.hdr           = hdr;

% filter data  
%------------
cfg            = [];
cfg.channel    = 'meg';         
cfg.detrend    = 'yes';
cfg.demean     = 'yes';        
cfg.continuous = 'yes';
cfg.bpfilter   = 'yes';
cfg.bpfreq     = bp_freq;       % e.g. [0.1 200];
data           = ft_preprocessing(cfg, data);
 
% check frequency-spectrum
%-------------------------
if check
    cfg             = [];
    cfg.method      = 'mtmfft';
    cfg.output      = 'pow';
    cfg.foi         = 0:1/2:200;
    cfg.taper       = 'hanning';
    cfg.pad         = 'maxperlen';
    freq_data       = ft_freqanalysis(cfg,data);

    channel2look = {'MT'};
    channel_out  = channel_selection(channel2look,freq_data.label);
    ind          = contains(freq_data.label,channel_out);
   
    figure('name','powerspectrum')
    plot(freq_data.freq,10*log10(freq_data.powspctrm(ind,:)));
    % xlim([0 50])
    title('filtered data')
    grid on
end

% downsample data
%----------------
cfg            = [];
cfg.resamplefs = fs_down;
cfg.detrend    = 'no';
data           = ft_resampledata(cfg,data);   
neurodata      = data.trial;


%% epoching and downsampling audiodata
%-------------------------------------
triallength       = info.triallength;
audiodata_epoched = cell(1,F);
downsampling      = 1;

for n = 1:F % files               
    audiodata_epoched{n} = epoch_audiodata(onset_envelope{n},fsamp,triallength,downsampling,fs_down);
end

% combine data
audiodata = [];
for n = 1:F
    audiodata = [audiodata,audiodata_epoched{n}.trials'];
end
% clear audiodata_epoched 

%% match length of audio and neuro epochs 
Nneuro = length(neurodata);
Naudio = length(audiodata);

if ~isequal(Nneuro,Naudio)
    error('error: number of epochs are different in neurodata and audiodata')
end

check_length = zeros(Naudio,2);
for n = 1:Naudio
    check_length(n,1) = size(neurodata{n},2);
    check_length(n,2) = length(audiodata{n});
end

% get index of minimum length for each trial
[M,~] = min(check_length,[],2);

% match length
for n = 1:Naudio
    neurodata{n} = neurodata{n}(:,1:M(n));
    audiodata{n} = audiodata{n}(1:M(n));
    data.time{n} = data.time{n}(1:M(n));
end

data.trial = neurodata;


%% save and segment data
% [audio_segments,neuro_segments,label_segments] = segment_data(audiodata,neurodata,winlen,shift,fs_down);
[audio_segments,neuro_segments,time_segments,tab_label_segments] = segment_data(audiodata,neurodata,data.time,winlen,shift,fs_down);
%%
crosscorr_segments.neurodata     = neuro_segments;
crosscorr_segments.audiodata     = audio_segments;

data.trial = neuro_segments;
data.time  = time_segments;

label = cat(1,data.label,'audio');
for n = 1:length(data.trial)
    triali{n} = cat(1,data.trial{n},audio_segments{n});
end
data.label = label;
data.trial = triali;
clear triali label

cfg         = [];
cfg.channel = 'audio';
audio1      = ft_selectdata(cfg,data);





info_previous      = info;
clear info

info.fsamp         = fs_down;
info.bp_freq_neuro = bp_freq;
info.lp_freq_audio = lp_freq_audio;
info.segmentlabel  = label_segments;
info.chanlabel     = data.label;
info.shift         = shift;
info.winlen        = winlen;
info.gradinfo      = data.grad;
info.info_previous = info_previous;

if save_data
    save(['helper_data' filesep filename2save '_Crosscorr_segments'],'Crosscorr_segments','-v7.3')
end



end

%% clean up
rmpath(fullfile('Z:','analysis','subject_files'))
rmpath(fullfile('Z:','analysis','preprocessing_batch','helper_functions'))
rmpath(fullfile('Z:','analysis','analysis_speech','helper_functions'));

%% functions
%--------------------------------------------------------------------------

function [envelope_diff] = cal_envelope(audiodata,fs,b,a)

% envelope
envelope = abs(hilbert(audiodata));
% Zero-phase digital filtering, lowpass 32 Hz
envelope_filtered = filtfilt(b,a,envelope);
% derivative (n-1) samples for onset envelopes
envelope_diff = diff(envelope_filtered)/(1/fs);
% add last sample artificially
envelope_diff = [envelope_diff; envelope_diff(end)];
% half-wave rectified
envelope_diff(envelope_diff < 0) = 0;

end

%--------------------------------------------------------------------------

function [audio_segments,neuro_segments,time_segments,tab_label_segments] = segment_data(audiodata,neurodata,timedata,winlen,shift,fsamp)
N_epochs = length(audiodata);

WindowLength = winlen*fsamp; % in samples
ShiftInt     = shift*fsamp;  % in samples

neuro_segments = [];
audio_segments = [];
time_segments  = [];
label_segments = zeros(N_epochs,2);
epoch_name     = [];
n              = 1;

for e = 1:N_epochs
    
    label_segments(e,1) = n;
    
    for s = 1:ShiftInt:(length(audiodata{e})-WindowLength+1)
      audio_segments{n} = audiodata{e}(s:s+WindowLength-1);
      neuro_segments{n} = neurodata{e}(:,s:s+WindowLength-1);
      time_segments{n}  = timedata{e}(s:s+WindowLength-1);
      n = n+1;
    end  
    
    add = {['epoch ',num2str(e)]};
    epoch_name = cat(1,epoch_name,add);
    
    label_segments(e,2) = n-1;
end

tab_label_segments = array2table(label_segments,...
                     'VariableNames',{'Start','End'},...
                     'Rownames',epoch_name);
end

%--------------------------------------------------------------------------

function scrollplot( dx, x, ax )
    % Set appropriate axis limits
    for ii = 1:numel(ax)
        set( ax(ii), 'xlim', [0 dx] );
    end

    % Create Uicontrol slider
    % The callback is another local function, this gives us more
    % flexibility than a character array.
    uicontrol('style','slider',...
        'units', 'normalized', 'position', [0.1 0.01 0.8 0.05],...
        'callback', @(slider, ~) scrollcallback( ax, dx, slider ), ...
        'min', 0, 'max', max(x)-dx );
end

%--------------------------------------------------------------------------

function scrollcallback( ax, dx, slider, varargin )
    % Scroller callback loops through the axes objects and updates the xlim
    val = slider.Value;
    for ii = 1:numel(ax)
        set( ax(ii), 'xlim', val + [0, dx] );
    end
end
