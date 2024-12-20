%--------------------------------------------------------------------------
% Till Habersetzer, 22.10.2024
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Computation of flat-spectrum chirps based on:
% - make_stim.m (Stefan Uppenkamp ~2000, creates flat spectrum chirps and
%   gates)
% - approximted chirp
%   Dau, Torsten, et al. "Auditory brainstem responses with optimized chirp
%   signals compensating basilar-membrane dispersion." The Journal of the 
%   Acoustical Society of America 107.3 (2000): 1530-1540.
% - phase-only reconstruction of the chirp signal (set amplitude to 1 for
%   all frequencies)
% 
%--------------------------------------------------------------------------

close all
clear all
clc

% creates flat spectrum chirps and gates
f1    = 100;	% f start
f2    = 10400;	% f stop
fsamp = 44100;	% sampling frequency
dur   = 0.012;	% duration of chirp (sec)
% dur   = 0;	% use minimum duration (sec)
alpha = 3.0;	% exponent alpha
lwin  = 0.001;	% duration of ramps (sec)

% call to bmchirp with var. alpha, find last non-zero sample

c_up = alchirp(f1,f2,fsamp,dur,alpha);
n = length(c_up);
for idx = n :-1:1
	if c_up(idx) ~= 0
		sample = idx;
		break
    end
end

% phase-only reconstruction of the chirp signal
fftc_up = fft(c_up);
fc_up   = real(ifft(cos(angle(fftc_up)) + 1i * sin(angle(fftc_up))));

% cut off zeros for gating, gate

c_up    = fc_up(1:sample);
n       = length(c_up);
nwin    = 2*round(fsamp*lwin);
window  = hanning(nwin,'symmetric')';
%window  = hamming(nwin)';
plateau = 1 + zeros(1,n-nwin);
gate    = [window(1:nwin/2), plateau, window(nwin/2+1:nwin)];
c_up    = c_up.* gate;
c_dn    = fliplr(c_up);

% show window
% window1 = hanning(nwin,'periodic')';
% window2 = hanning(nwin,'symmetric')';
% figure
% plot(window1)
% hold on
% plot(window2)
% legend('periodic','symmetric')

% write wav files
% normalize
if ~exist('final_stimuli','dir')
    mkdir final_stimuli
end
audiowrite(fullfile('final_stimuli','up.wav'),c_up/max(abs(c_up)),fsamp);
audiowrite(fullfile('final_stimuli','down.wav'),c_dn/max(abs(c_up)),fsamp);

% soundsc(c_up,fsamp)

% figure
% plot(c_up)

%% Add click
%--------------------------------------------------------------------------
dur = 100*10^-6; % 100 us
% shorten click
len = ceil(dur*fsamp);

click      = ones(1,len);
click      = [0,click,0];
% sound(click,fsamp)
audiowrite(fullfile('final_stimuli','click.wav'),click,fsamp);

% [click,fsamp] = audioread('click.wav');
% audiowrite(['final_stimuli',filesep,'click.wav'],click(1:4)/max(abs(click)),fsamp);

