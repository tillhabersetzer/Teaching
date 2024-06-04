close all;
clear all;
clc;

% The stimuli are rescaled so that they all have the same rms value within
% the experimental paradigm given a specific repitition rate.
p = fileparts(mfilename('fullpath'));

% Settings for Calibration Stimuli
%--------------------------------------------------------------------------

% Import Settings
Settings = importdata(fullfile('..','Settings','Settings.mat'));

% Generate calibration signal, concatenate 180 Up-chirps according to 
% stimulus paradgim (3 chirps per sec x 10 ~ 10s)
Nstim  = 30;
Jitter = Settings.Jitter;
% Set the random number generator to the default seed (0) 
rng('default')
s = rng; % save generator settings;
save(fullfile(p,'generator_settings.mat'))
Jitter = Jitter(1) + Jitter(2)*rand(1,Nstim) ; % [0.35,0.4]s

% Generate three Calibration Signals (Click, Up-Chirp, Down-Chirp) to
% calculate RMS value
%--------------------------------------------------------------------------

[click,SampFreq] = audioread(fullfile(p,'..','Stimuli','raw','click.wav'));% 1
[up,~]           = audioread(fullfile(p,'..','Stimuli','raw','up.wav'));   % 2
[down,~]         = audioread(fullfile(p,'..','Stimuli','raw','down.wav')); % 3
Jitter           = round(Jitter*SampFreq); % in samples

cal_signal = cell(1,3);
for n=1:Nstim
    cal_signal{1} = vertcat(cal_signal{1},vertcat(click,zeros(Jitter(n),1)));
    cal_signal{2} = vertcat(cal_signal{2},vertcat(up,zeros(Jitter(n),1)));
    cal_signal{3} = vertcat(cal_signal{3},vertcat(down,zeros(Jitter(n),1)));
end

% duration in sec
%dur = length(cal_signal)/SampFreq;

% Scaling: Click is reference due to lowest SNR (can't be upscaled)
% rms_click = scaling * rms_up -> scaling = rms_click/rms_up / rms_up = rms_down
rms_click = rms(cal_signal{1});
rms_up    = rms(cal_signal{2});
%rms_down  = rms(cal_signal{3});
scaling   = rms_click/rms_up;

% Rescale Signals and check RMS
rms_click_rs = 20*log10(rms(cal_signal{1}));
rms_up_rs    = 20*log10(rms(scaling*cal_signal{2}));
rms_down_rs  = 20*log10(rms(scaling*cal_signal{3}));

% Save Signals
%--------------------------------------------------------------------------
audiowrite(fullfile('..','Stimuli','rms_equalized','click_rmseq.wav'),click,SampFreq);
audiowrite(fullfile('..','Stimuli','rms_equalized','up_rmseq.wav'),scaling*up,SampFreq);
audiowrite(fullfile('..','Stimuli','rms_equalized','down_rmseq.wav'),scaling*down,SampFreq);

% save rescaled up-chirp as calibration signal
audiowrite(fullfile('..','Stimuli','calibration','calibration_up_rmseq.wav'),scaling*cal_signal{2},SampFreq);

% figure
%plot(scaling*cal_signal{2})

% listen
% sound(cal_signal{1},SampFreq)