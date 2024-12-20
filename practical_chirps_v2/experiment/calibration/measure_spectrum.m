%--------------------------------------------------------------------------
% @author: Till Habersetzer
%          Carl von Ossietzky University Oldenburg
%          till.habersetzer@uol.de
%
%          Created: 29.10.24
%          Finished first draft: 02.11.24
%          
%          Revisions: 
%
% This script estimates and plots the power spectrum for a sequence of
% transient signals.
% A sequence of transient stimuli is played back, and the microphone signal
% is recorded. Based on this recording, a spectrum is estimated for each
% transient pulse and for the entire sequence by averaging across all 
% epochs.
%--------------------------------------------------------------------------

close all
clearvars
clc

%% 1.) Playback and record transient
%--------------------------------------------------------------------------
% Audio output channel 1 is mapped to audio input channel 1 to create a 
% trigger signal, which is recorded as a loopback signal. This serves as a 
% reference for the transient sequence.
filename                  = 'click_eq.wav';
channel_out               = [0,1]; % audio, trigger
channel_in                = [0,1]; % microphone, trigger loopback
gain_dB                   = -20;
n_stim                    = 100;
gap                       = 1; % 1s
show_audio_tracks         = true;
[rec_filenames,trial_dur] = playrec_stim(channel_out, channel_in, gain_dB, filename, n_stim, gap, show_audio_tracks);

%% 2.) Compute spectrum
%--------------------------------------------------------------------------

% Load and plot signals
%----------------------
% rec_filenames = {'up_recording.wav','up_trigger.wav'};
[mic_sig, fs]    = audioread(rec_filenames{1});
[trigger_sig, ~] = audioread(rec_filenames{2});

figure
plot(mic_sig./max(abs(mic_sig)))
hold on
plot(trigger_sig./max(abs(trigger_sig)))
legend('Microphone','Trigger')

%% 2.1) Estimate latency between microphone signal and transient (click) sequence
%--------------------------------------------------------------------------
% based on peaks in recorded signals.

% Detect peaks in channels and plot them
%---------------------------------------
% normalized derivative
deriv_trigger = diff(trigger_sig)/max(abs(diff(trigger_sig)));
deriv_mic     = diff(mic_sig)/max(abs(diff(mic_sig)));

[~, peak_idx_trigger] = findpeaks(deriv_trigger, 'MinPeakHeight', 0.5, 'MinPeakDistance', trial_dur-0.1*fs, 'NPeaks', n_stim);
[~, peak_idx_mic]     = findpeaks(deriv_mic, 'MinPeakHeight', 0.5, 'MinPeakDistance', trial_dur-0.1*fs, 'NPeaks', n_stim);

if length(peak_idx_trigger) ~= n_stim || length(peak_idx_mic) ~= n_stim
    warning('Not all peaks detected! / Potential mismatch in detected peaks!')
end

% Visualize
%----------
timy = 1:length(trigger_sig);
figure
subplot(2,1,1)
plot(timy,trigger_sig)
hold on
plot(timy(peak_idx_trigger),trigger_sig(peak_idx_trigger),'ro')
title('Trigger signal')
subplot(2,1,2)
plot(timy,mic_sig)
hold on
plot(timy(peak_idx_mic),mic_sig(peak_idx_mic),'ro')
title('Microphone signal')

% Estimate latency based on peaks
%--------------------------------
peaks_latency     = peak_idx_mic - peak_idx_trigger;
peaks_latency_avg = ceil(mean(peaks_latency));

figure
histogram(peaks_latency/fs)
hold on
yLimits = ylim;
plot([peaks_latency_avg peaks_latency_avg]/fs, yLimits, 'r--', 'LineWidth', 2);
title('Latency / s')

%% 2.2) Estimate spectrum
%--------------------------------------------------------------------------

% Cut into epochs and compute spectrum
%-------------------------------------
% 1. either based on detected peaks in microphone signal
peak_idx = peak_idx_mic;
% 2. or based on trigger onsets with latency correction
% peak_idx = peak_idx_trigger + peaks_latency_avg;

% Epoch window is based on gap size
lower_bound = ceil(0.01*trial_dur); % 1 % of gap
upper_bound = ceil(0.1*trial_dur); % 10 % of gap

n_fft    = 2^nextpow2(lower_bound+upper_bound+1); % zeropad to 1s, frequency resolution 1Hz
L        = n_fft;
freqs    = fs*(0:(L/2))/L;
spectrum = zeros(length(freqs),n_stim);

for nidx = 1:n_stim
    % Epoching
    epoch = mic_sig(peak_idx(nidx)-lower_bound:peak_idx(nidx)+upper_bound);

    % Spectrum
    L                = n_fft;
    Y                = fft(epoch,n_fft);
    P                = abs(Y/L).^2; % normalize and power
    P                = P(1:L/2+1);
    spectrum(:,nidx) = P;
    clear epoch L Y P

end

% Compute mean spectrum
spectrum_avg = mean(spectrum,2);

% Plot spectrum
%--------------
figure
for nidx = 1:n_stim
    plot(freqs/1000,db(spectrum(:,nidx),'power')) % dB: 10*log10(x)
    hold on
end
% Add mean
plot(freqs/1000,db(spectrum_avg,'power'),'k','linewidth',2) % dB: 10*log10(x)
xlabel("f / kHz")
ylabel("|P1(f)|^2 / dB")
xlim([0.01,20])
% xlim([f1*0.001/10,2*f2*0.001])
title(sprintf('Spectrum | epochs: %i', n_stim),'fontsize',16)
grid('on')
axs = gca;
set(axs,'xtick', [10,100,1000,5000,10000,15000]/1000, 'xticklabel',{'0.01' '0.1' '1' '5' '10' '15'},'fontsize',fontsize,'LineWidth',linewidth)
set(axs,'XScale', 'log');

%% Functions
%--------------------------------------------------------------------------
function [rec_filenames,trial_dur] = playrec_stim(channel_out, channel_in, gain_dB, filename, n_stim, gap, show_audio_tracks)
%--------------------------------------------------------------------------
% This function plays a calibration stimulus while simultaneously recording 
% the output. The stimulus, which is a transient calibration signal, is 
% looped a specified number of times. This way, a sequence is generated.
% The gap size between consecutive elements and its number can be defined.
% It uses the SoundMexPro library for audio playback and recording, and 
% optionally visualizes audio tracks.
% 
% Additionally to the calibration signals, a trigger signal is played
% back on the sedond audio channel which is aligned to the calibration
% signal.
%
% Wiring:
% Connect the second output channel to the second input channel for the
% trigger (loopback).
% 
% Parameters:
%------------
% channel_out (array of integers): Output channel for playback on the sound card.
%   0: AnalogOut 0 - audio signal
%   1: AnalogOut 1 - trigger signal
%   Example: [0, 1] - where AnalogOut 0 is a microphone signal, and AnalogOut 1 is a loopback trigger signal.
%
% channel_in (array of integers): Input channels for recording from the sound card.
%   0: AnalogIn 0 - audio signal
%   1: AnalogIn 1 - trigger signal
%   Example: [0, 1] - where AnalogIn 0 is a microphone signal, and AnalogIn 1 is a loopback trigger signal.
%
% gain_dB (float): Gain applied to the audio signal in decibels (dB).
%   Positive values increase gain (amplification).
%   Negative values decrease gain (attenuation).
%
% filename (string): Filename of the calibration signal file to be played back 
%   (e.g., up.wav).
%
% n_stim (integer, optional): Number of times to loop the calibration signal.
% 
% gap (float, optional): Duration of the gap (in seconds) between successive stimuli.
%
% show_audio_tracks (boolean, optional): If true, visualizes the playback and recording tracks using SoundMexPro.
%
% Returns:
%---------
% rec_filenames (cell array of strings): 
%   Names of the recorded audio files (one for each input channel).
% trial_dur (integer): 
%   Duration of a single stimulus trial, calculated as the length of the 
%   audio signal plus the gap (in samples).
%--------------------------------------------------------------------------

% Set default function parameters
%--------------------------------
if nargin<7
    show_audio_tracks = false;
end
if nargin<6
    gap = 1; % 1 s
end
if nargin<6
    n_stim = 10; % 10 repititions
end
if nargin<4
    filename = 'up.wav';
end
if nargin<3
    gain_dB = 0; % 0 gain
end
if nargin<2
    channel_in = [0,1]; % second channel used for trigger
end
if nargin<2
    channel_out = [0,1] ; % second channel used for trigger
end

fprintf('\nMeasure spectrum of transients\n')
fprintf('------------------------------\n\n')

% Settings
%--------------------------------------------------------------------------
% Import Settings
settings = importdata(fullfile('..','settings','settings_hearing_threshold.mat'));

% Import audio file
[audio, fs] = audioread(fullfile('..','audio','tip300',filename));
% [audio, fs] = audioread(fullfile('..','audio','sensimetrics',filename));
if settings.samp_freq~=fs
    error("Sampling frequencies don't match (%i~=%i).",settings.samp_freq,fs)
end

% Extract correct channel if stereo
% if size(audio,2)>1
%     audio = audio(:, 2);
% end

% Get taskname for file naming
[~, taskname, ~] = fileparts(filename);

% Generate signal
%----------------
gap_s     = ceil(fs*gap); % gap between stimuli in samples
trial_dur = length(audio)+gap_s; % return duration of single trial
polarity  = 1;
gap_s     = ceil(fs*gap); % gap between stimuli in samples
% Polarity of audio is inversed analogously to experiment

trig                         = zeros(size(audio));
trig(1:ceil(length(trig)/2)) = 0.4;

% Add 1s zeros
cal_audio = vertcat(zeros(fs,1),repmat(vertcat(polarity*audio,zeros(gap_s,1)),n_stim,1));
trigger   = vertcat(zeros(fs,1),repmat(vertcat(trig,zeros(gap_s,1)),n_stim,1));

% figure
% hold on
% plot(cal_audio)
% plot(trigger)
% legend('audio','trigger')

% SoundmexPro
%--------------------------------------------------------------------------
% Attenuated calibration signal
% level_stim_new = gain_dB + level_stim
% -> scaling = 10^(gain_dB/20)
cal_audio = 10.^(gain_dB/20)*cal_audio;

% check for clipping
if max(abs(cal_audio))>1
    error('Signal is clipping!')
end

% Level of initial calibration signal in dB FS (peak)
level_stim_init = 20*log10(max(abs(audio)));
% Level of scaled signal in dB FS (peak)
level_stim = 20*log10(max(abs(cal_audio)));

fprintf('Initial Stimulus Level: %.1f dB FS (peak).\n',level_stim_init)
fprintf('Stimulus Level is changed by %.1f dB.\n', gain_dB)
fprintf('Adjusted stimulus Level: %.1f dB FS (peak).\n', level_stim)

ok = soundmexpro('init', ...            % command name
    'force',      1, ...                % exit internally called before init
    'driver',     settings.asio_driver, ...
    'samplerate', settings.samp_freq, ...
    'output',     channel_out, ... % add additional trigger channel
    'input',      channel_in, ...
    'track',      2, ...  
    'autocleardata', 1);            % clear audio data on memory
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

% set filename of recording of first channel to other than default and
% print both resulting filenames
[ok, rec_filenames] = soundmexpro('recfilename', ...   
    'filename', {sprintf('%s_recording.wav',taskname), sprintf('%s_trigger.wav',taskname)}, ...               
    'channel', channel_in ...                                % channel to set filename
    );
if ~ok, error(['error calling ''recfilename''' error_loc(dbstack)]); end
fprintf('First channel does record to file: %s\n', rec_filenames{1})
fprintf('Second channel does record to file: %s\n', rec_filenames{2})

ok = soundmexpro('loadmem', ...     % command name
    'data', [cal_audio,trigger], ...   % data vector
    'loopcount', 1);                % 1: only once
if ~ok, error('error calling ''loadmem'''); end

% show visualization
if show_audio_tracks
    ok =  soundmexpro('showtracks');
    if ~ok,error(['error calling ''showtracks''' error_loc(dbstack)]); end
end

% Playback and recording!
ok = soundmexpro('start', ...
    'length',0); % device is never stopped, zeros are played endlessly
if ~ok,  error('error calling ''start''' ); end

% Now show, that both channels are recording at the moment
[ok, recording] = soundmexpro('recording');
if ok ~= 1
    error(['error calling ''recording''' error_loc(dbstack)]);
end
fprintf('First channel does record: %i\n', recording(1))
fprintf('Second channel does record: %i\n', recording(2))

playing = 1;
while max(playing)
    % query device
    [ok, playing] = soundmexpro('playing');
    % check success of command itself
    if ok ~= 1
        clear soundmexpro
        error(['error calling ''playing''' error_loc(dbstack)]);
    end
    % NOTE: no while loops without a small pause!
    pause(0.01);
end

% Add additiobal second
pause(1)

% show identical recpositions (recorded lengths) of the two channels
[success, recpos] = soundmexpro('recposition');
if success ~= 1
    error(['error calling ''recpause''' error_loc(dbstack)]);
end
fprintf('Samples recorded on first channel: %i\n',recpos(1))
fprintf('Samples recorded on second channel: %i\n',recpos(2))

% Check for Clipping
[success, clipout, clipin] = soundmexpro('clipcount');
if success ~= 1
    error(['error calling ''clipcount''' error_loc(dbstack)]);
end
fprintf('Clipped input buffers on first channel: %i.\n',clipin(1)) 
fprintf('Clipped input buffers on second channel: %i.\n',clipin(2))
fprintf('Clipped output buffers on first channel: %i.\n',clipout(1)) 
fprintf('Clipped output buffers on second channel: %i.\n',clipout(2))

ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end 

end
