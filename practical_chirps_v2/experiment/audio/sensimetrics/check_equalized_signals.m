close all
clearvars
clc

condition = 'click';
% condition = 'up';
% condition = 'down';
% condition = 'sinus_1000Hz';

% Import left right and both
[audio_left,fs]  = audioread([condition,'_eq_left.wav']);
[audio_right,fs] = audioread([condition,'_eq_right.wav']);
[audio_both,fs]  = audioread([condition,'_eq_both.wav']);

% compare signals
if ~isequal(audio_left(:,1),audio_left(:,2))
    warning("Left earphone: Left and right signals aren't the same!")
end
if ~isequal(audio_right(:,1),audio_right(:,2))
    warning("Right earphone: Left and right signals aren't the same!")
end
if isequal(audio_both(:,1),audio_both(:,2))
    warning("Both earphones: Left and right signals are the same!")
end

time = (0:length(audio_left)-1)/fs;

figure
subplot(2,1,1)
plot(time,audio_left(:,1),'x-')
hold on
plot(time,audio_both(:,1),'x-');
legend('left channel only', 'both channels - left')
title('Left')
subplot(2,1,2)
plot(time,audio_right(:,2),'x-')
hold on
plot(time,audio_both(:,2),'x-');
legend('right channel only', 'both channels - right')
title('Right')