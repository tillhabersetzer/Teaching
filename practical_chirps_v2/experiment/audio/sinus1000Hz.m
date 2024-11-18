%--------------------------------------------------------------------------
% @author: Till Habersetzer
%          Carl von Ossietzky University Oldenburg
%          till.habersetzer@uol.de
%          
%          Created: 27.10.24
%--------------------------------------------------------------------------

close all
clear all
clc

% Import settings
settings = importdata(fullfile('..','settings','settings_hearing_threshold.mat'));
fs       = settings.samp_freq;

% Generate 1000 Hz sinus (1 sec) and concatenate
freq = 1000;
% time = 0:1/settings.samp_freq:1-1/settings.samp_freq;
% time = (0:fs-1)/fs; % 1s
time = (0:10*fs-1)/fs; % 10s
% time = (0:fs/freq-1)/fs; % only single period
sig  = sin(2*pi*freq*time);

% Save audio
%-----------
audiowrite('sinus_1000Hz.wav',sig,fs);

% Check that transition works
%----------------------------
figure
plot(time,sig,'-x')
hold on
plot(time+time(end)+1/fs,sig,'-x')
xlim([1-1/freq,1+1/freq])

% Soundcheck
%-----------
audio = [];
for n=1:2
    audio = vertcat(audio,sig');
end

soundsc(audio,fs)

