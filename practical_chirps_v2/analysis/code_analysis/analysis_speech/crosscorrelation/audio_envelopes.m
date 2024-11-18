close all; clear all; clc

% choose envelope type for crosscorrelation
envelopetype{1} = 'onset_envelope';
envelopetype{2} = 'envelope';

% lowpass audiodata
lp_freq_audio = 25;

% addpath for basic functions
addpath(fullfile('Z:','analysis','analysis_speech','helper_functions'));

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

% hilbert/filter/derivative/half wave rectified
%----------------------------------------------
envelopes       = cell(2,4);
raw_audiodata   = cell(1,4);
for n = 1:4
    [raw_audiodata{n},fs] = audioread(path_audiodata{n}); 
    if ~isequal(fs,fsamp)
        error('unexpected sampling frequency!')
    end
    clear fs
    for t = 1:2
        envelopes{t,n} = cal_envelope(raw_audiodata{n},fsamp,b,a,envelopetype{t});
        % zscore data
%         envelopes{t,n} = zscore(envelopes{t,n},0,1);
    end
end



% min(envelopes{1,4})

%% have a look at speech onset envelope
%--------------------------------------------------------------------------
% choose envelope type to plot
% 1.) onset envelope 2.) "normal envelope"
type = 2;

label = {'Das_schwatzende_Herz1','Das_schwatzende_Herz2',...
         'Die_Maske_des_Roten_Todes1','Die_Maske_des_Roten_Todes2'};
% choose file (n=1,2,3)
n = 2;
x = (1:length(raw_audiodata{n}))*(1/fsamp);
% dx is the width of the axis 'widxow'
dx = 20; % in sec
% Initialise the figure once, and we only need to set the properties once
fig{t} = figure; clf;
set( fig{t}, 'doublebuffer', 'on'); 
% Create a placeholder for axes objects
ax = gobjects(2,1);
% Create plots, storing them in the axes object
ax(1) = subplot(3,1,1);
plot(x, raw_audiodata{n});
title('raw audiodata')
ax(2) = subplot(3,1,2);
plot(x,envelopes{type,n}); 
xlabel('t / s')
title('speech onset envelope')
ax(3) = subplot(3,1,3);
plot(x, raw_audiodata{n});
hold on
plot(x,envelopes{type,n},'r'); 
xlabel('t / s')
title('together')
% Set up the scroller for the array of axes objects in 'ax'
scrollplot( dx, x, ax)
sgtitle([envelopetype{type},' : ',label{n}]) 

%--------------------------------------------------------------------------

%% Clean up 
rmpath(fullfile('Z:','analysis','analysis_speech','helper_functions'));

%% functions
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

%--------------------------------------------------------------------------


