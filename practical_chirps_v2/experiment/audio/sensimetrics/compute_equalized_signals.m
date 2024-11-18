close all
clearvars
clc

% Add filter function
root_path = fullfile('..','..','..','..','transfer_function','sensimetrics_s14');
addpath(fullfile(root_path,'EQfiltering_Matlab_Utility'))
% addpath('C:\Users\tillhabersetzer\Nextcloud\Synchronisation\Projekte\GitHub\Teaching\practical_chirps_v2\transfer_function\sensimetrics_s14\EQfiltering_Matlab_Utility')

% Load filter
%------------
% left
[h_left,fs1] = load_filter(fullfile(root_path,'600L_600R','EQF_600L.bin'));
% right
[h_right,fs2] = load_filter(fullfile(root_path,'600L_600R','EQF_600R.bin'));

if ~isequal(fs1,fs2)
    error('Sampling frequencies are different!')
else
    fs = fs1;
end

% reduce size in samples
% lower_bound  = 214;
% upper_bound = 318;
% h_left  = h_left(lower_bound:upper_bound);
% h_right = h_right(lower_bound:upper_bound);

time = (0:length(h_left)-1)/fs;
figure
subplot(2,2,1)
plot(time,h_left,'bx-')
title('Left impulse response')
subplot(2,2,2)
plot(time,h_right,'rx-')
title('Right impulse response')
subplot(2,2,3:4)
plot(h_left,'bx-')
hold on
plot(h_right,'rx-')
title('Both together')
legend('Left impulse response','Right impulse response')

%% Load audio 
%--------------------------------------------------------------------------
condition = 'click';
% condition = 'up';
% condition = 'down';
% condition = 'sinus_1000Hz';

% Import left right and both
[audio,fs]  = audioread(fullfile('..','tip300',[condition,'.wav']));

if ~isequal(fs,fs1)
    error('Sampling frequencies are different!')
end

if size(audio,2)>1
    error('No mono signal!')
end

figure
plot(audio,'x-')

%% Apply filter
%--------------------------------------------------------------------------
% audio_left_eq  = conv(audio, h_left, 'same');  % Apply filtering with 'same' length
% audio_right_eq = conv(audio, h_right, 'same');  % Apply filtering with 'same' length

audio_left_eq  = conv(audio, h_left);  
audio_right_eq = conv(audio, h_right);  

figure
subplot(3,2,1)
plot(audio,'x-')
title('Audio')
subplot(3,2,2)
plot(time,h_left,'bx-')
hold on
plot(time,h_right,'rx-')
title('Both together')
legend('Left impulse response','Right impulse response')
subplot(3,2,3)
plot(audio_left_eq,'x-')
title('Audio left equalized')
subplot(3,2,4)
plot(audio_right_eq,'x-')
title('Audio right equalized')
subplot(3,2,5:6)
plot(audio_left_eq,'bx-')
hold on
plot(audio_right_eq,'rx-')
legend('Audio left equalized','Audio right equalized')
grid on

%% Cut signals and multiply with ramp to correct for onset effects....
%--------------------------------------------------------------------------

% lower_bound  = 214;
% upper_bound = 318;

% Cut signal
%-----------
audio_left_eq_cut  = audio_left_eq;
audio_right_eq_cut = audio_right_eq;

% Apply gate
%-----------
lwin              = 0.001;
n                 = length(audio_left_eq_cut);
nwin              = 2*round(fs*lwin);
window            = hanning(nwin,'symmetric')';
plateau           = 1 + zeros(1,n-nwin);
gate              = [window(1:nwin/2), plateau, window(nwin/2+1:nwin)];
audio_left_final  = audio_left_eq_cut .* gate';
audio_right_final = audio_right_eq_cut .* gate';

% Show signals
%-------------
figure
subplot(2,2,1)
plot(audio_left_final,'x-')
title('Audio left equalized')
subplot(2,2,2)
plot(audio_right_final,'x-')
title('Audio right equalized')
subplot(2,2,3:4)
plot(audio_left_final,'bx-')
hold on
plot(audio_right_final,'rx-')
legend('Audio left equalized','Audio right equalized')

%% Save audio file
%--------------------------------------------------------------------------
audio = [audio_left_final./max(abs(audio_left_final)),audio_right_final./max(abs(audio_right_final))];

audiowrite([condition,'_eq.wav'],audio,fs);










