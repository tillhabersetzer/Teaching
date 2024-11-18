%--------------------------------------------------------------------------
% @author: Till Habersetzer
%          Carl von Ossietzky University Oldenburg
%          till.habersetzer@uol.de
%          
%          Revisions: 
%          - 16.07.2024: Implementation in Python
%          - 25.10.2024: Comversion to Matlab
% 
% Calibration
%--------------------------------------------------------------------------
%
% Best practice:
%---------------
% - Playback Transient: 
%   Play back a transient signal with a e.g. 30 dB software attenuation
% - Hardware Gain Adjustment: 
%   Adjust the hardware gain so that the sound level meter (SLM) reads 80 dB.
% - Maximum Possible Level: 
%   This setup allows for a maximum possible level of 110 dB without attenuation.
% 
% Level Definitions:
%-------------------
% - level_transient: 
%   Level of the transient signal in dB FS, measured as the peak value
%   max(abs(...))
% - level_slm: 
%   Level measured by the sound level meter, optimally as dB (p-p) peSPL.
%   Adjust 1000 Hz sine wave with a peak-to-peak voltage (ppv) that is 
%   equal to the ppv of the transient and measure its exponential 
%   time-weighted sound level 
%
%   see: 
%   Laukli, Einar, and Robert Burkard. "Calibration/standardization of 
%   short-duration stimuli." Seminars in hearing. Vol. 36. No. 01. 
%   Thieme Medical Publishers, 2015.
%   
% Calibration Value Calculation:
%-------------------------------
% 
% - Without Attenuation:
%   cal_val = level_slm - level_transient
% 
% - With Attenuation (gain_dB):
%   cal_val = level_slm - (level_transient + gain_dB)
% 
% Interpretation:
%----------------
% The cal_val represents the sound level meter reading when the signal is 
% played at full amplification (0 dB FS). The maximum achievable level 
% should be sufficiently high.
% 
% Example:
% level_transient = -20 dB FS
% level_slm       = 80 dB
% cal_val         = 100 dB 
%
% This means that if the transient signal were played at 0 dB FS, the sound 
% level meter would read 100 dB.
% 
% Recommended Calibration Equipment for Insert Earphones
% ------------------------------------------------------
% - B&K Type 4157 Ear Simulator
% - B&K External Ear Simulator DB2012 with secure cap
% - B&K 2250 Sound Level Meter/Analyzer
% - B&K ZC-0032 microphone preamplifier 
% - B&K Type 4231 Sound Calibrator  
%--------------------------------------------------------------------------

close all
clearvars
clc

%% Calibration
%--------------------------------------------------------------------------

%% 1.) Playback transient
% - Measure p-p voltage on oscilloscope
% - Note down stim_dB_Fs_transient
filename             = 'up';
channel              = 1;
gain_dB              = -20;
show_audio_tracks    = false;
level_stim_transient = calibration_transient(channel, gain_dB, filename, show_audio_tracks);

%% 2.) Playback sinus
% - Adjust the gain_dB so that the p-p voltage of the sinus matches the
%   p-p voltage of the transient
% - Note down the level on the sound level meter (level_slm)
filename          = 'sinus';
channel           = 1;
gain_dB           = -15;
show_audio_tracks = false;
level_stim_sinus  = calibration_transient(channel, gain_dB, filename, show_audio_tracks);

%% 3.) Compute the calibration value
% cal_val = level_slm - level_stim_transient

%% Function
%--------------------------------------------------------------------------
function [level_stim] = calibration_transient(channel, gain_dB, filename, show_audio_tracks)
%--------------------------------------------------------------------------
% Playback of calibration signal. The calibration signal (e.g. transient)
% is concatenated in an infinite loop. 
% Therefore, the transient is stretched with 0s so that it has 1 sec duration.
% Transient: Calibration with max-abs-value
% 
% Parameters:
%------------
% channel: integer
%   Output channel of sound card for playback.
%   0: AnalogOut 0 - left 
%   1: AnalogOut 1 - right
% 
% gain_dB: float
%   Channel gain in dB.
%   >0: amplification
%   <0: attenuation
%
% filename: str
%   Filename of calibration signal
%   'up'   : up-chirp
%   'click': click
%   'dowm' : down-chirp
%   'sinus': 1000 Hz sinus
%
% show_audio_tracks: bool
%   Shows visualization of files/vectors in tracks with SoundMexPro.
%   True
%   False
% 
% Returns:
%---------
% level : float in dB FS
% 
%--------------------------------------------------------------------------

% Set default function parameters
%--------------------------------
if nargin<4
    show_audio_tracks = false;
end
if nargin<3
    filename = 'up';
end
if nargin<2
    gain_dB = 0; % 0 gain
end
if nargin<1
    channel = 0; % left channel
end

fprintf('\nCalibration of Transient\n')
fprintf('------------------------\n\n')
fprintf('\nPlease not down calibration equipment and calibration values!\n')
fprintf('Press ESCAPE to stop playback.\n\n') 

% Settings
%--------------------------------------------------------------------------
% Import Settings
settings = importdata(fullfile('..','settings','settings_experiment.mat'));

% Import calibration signal
switch filename
    case 'up'
        audiofile = 'up.wav';
        % audiofile = 'up_eq.wav';
    case 'down'
        audiofile = 'down.wav';
        % audiofile = 'down_eq.wav';
    case 'click'
        audiofile = 'click.wav';
        % audiofile = 'click_eq.wav';
    case 'sinus'
        audiofile = 'sinus_1000Hz.wav'; % 10s duration
        % audiofile = 'sinus_1000Hz_eq.wav';
    otherwise
        error('Unexpected audio filename: %s!',filename)
end

[audio, fs] = audioread(fullfile('..','audio','tip300',audiofile));
% [audio, fs] = audioread(fullfile('..','audio','sensimetrics',audiofile));

% Extract correct channel if stereo
if size(audio,2)>1
    audio = audio(:, 2);
end

if settings.samp_freq~=fs
    error("Sampling frequencies don't match (%i~=%i).",settings.samp_freq,fs)
end

% Append click, up, down with zeros to 1s
%----------------------------------------
switch filename
    case {'up','down','click'}
        audio = vertcat(audio,zeros(fs-length(audio),1));
end

% SoundmexPro
%--------------------------------------------------------------------------
% Attenuated calibration signal
% level_stim_new = gain_dB + level_stim
% -> scaling = 10^(gain_dB/20)
calsig = 10.^(gain_dB/20)*audio;

% Level of initial calibration signal in dB FS (peak)
level_stim_init = 20*log10(max(abs(audio)));
% Level of scaled signal in dB FS (peak)
level_stim = 20*log10(max(abs(calsig)));

fprintf('Initial Stimulus Level: %.1f dB FS (peak).\n',level_stim_init)
fprintf('Stimulus Level is changed by %.1f dB.\n', gain_dB)
fprintf('Adjusted stimulus Level: %.1f dB FS (peak).\n', level_stim)

ok = soundmexpro('init', ...            % command name
    'force',      1, ...                % exit internally called before init
    'driver',     settings.asio_driver, ...
    'samplerate', settings.samp_freq, ...
    'output',     channel, ...
    'input',      [], ...
    'autocleardata', 1);            % clear audio data on memory
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

ok = soundmexpro('loadmem', ...     % command name
    'data', calsig, ...             % data vector
    'loopcount', 0);                % 0: endless loop 
if ~ok, error('error calling ''loadmem'''); end

% show visualization
if show_audio_tracks
    ok =  soundmexpro('showtracks');
    if ~ok,error(['error calling ''showtracks''' error_loc(dbstack)]); end
end

ok = soundmexpro('start', ...
    'length',0); % device is never stopped, zeros are played endlessly
if ~ok,  error('error calling ''start''' ); end

answer = questdlg("Press 'Yes' to stop the playback.",'Input','Yes','Yes');
switch answer
    case 'Yes'
        ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
        ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end 
end

end
