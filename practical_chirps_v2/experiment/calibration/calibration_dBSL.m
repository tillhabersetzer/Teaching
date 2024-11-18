function [threshold_estimate] = calibration_dBSL(subjectname,channel,filename)
%--------------------------------------------------------------------------
% % Till Habersetzer, 28.10.2024
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de
%
% CALIBRATION_DBSL Measures hearing threshold adaptively using transient 
% signals (e.g. up-chirps, clicks).
% This function presents a series of transients with a specified inter-stimulus 
% interval (ISI) and jitter, adjusting the gain until the hearing threshold 
% is determined. The signal can be presented to one ear (monaurally) or 
% both ears (diotically).
%
% The initial step size for gain adjustment is 10 dB, which decreases to 
% 5 and 2 dB as the procedure progresses through reversal points.
%
% Output:
%   - threshold_estimate: Adjusted signal level in in (p-p) peSPL that 
%                         corresponds to the determined hearing threshold.
% 
% Input:
%   - subjectname: Name of the subject (e.g., 'sub-01').
%   - channel: 
%       0 - left ear,
%       1 - right ear,
%       [0, 1] - both ears.
%   - filename: Name of the audio file ('click', 'up', or 'down') used for 
%               calibration.
%
% Example:
% threshold_estimate = calibration_dBSL('sub-01', 0, 'up');
%--------------------------------------------------------------------------

% TO DO
% add limitation of gain: reasonable boundaries

%% Settings
%--------------------------------------------------------------------------
close all % close other figures

% Set default function paramters
if nargin<3
    filename = 'up';
end
if nargin<2
    channel = 0;
end
if nargin<1
    subjectname = 'JohnDoe';
end

% Import Settings
settings = importdata(fullfile('..','settings','settings_hearing_threshold.mat'));
% Add avoidOverwrite function
addpath(fullfile('..','additional_functions')); 

% Import calibration signal
switch filename
    case 'up'
        audiofile = 'up.wav';
        cal_val   = settings.calibration.cal_val_up;
        start_val = settings.calibration.start_val_up;
    case 'down'
        audiofile = 'down.wav';
    case 'click'
        audiofile = 'click.wav';
        cal_val   = settings.calibration.cal_val_click;
        start_val = settings.calibration.start_val_click;
    otherwise
        error('Unexpected audio filename: %s!',filename)
end

[audio, fs] = audioread(fullfile('..','audio',audiofile));
% audio = audio./10;
if settings.samp_freq~=fs
    error("Sampling frequencies don't match (%i~=%i).",settings.samp_freq,fs)
end

% Generate calibration signal, concatenate ~180 stimuli according 
% to stimulus paradgim (3 chirps per sec x 10 ~ 10s)
%--------------------------------------------------------------------------
n_stim = 30; % must be an even number if multiple blocks are concatenated with reveres polarity
jitter = settings.jitter;

% Define seed, always the same randomization
rng(13)
jitter = jitter(1) + jitter(2)*rand(1,n_stim) ; % [0.35,0.4] s
jitter = round(jitter*fs); % in samples

cal_audio = [];
polarity  = 1;
% Polarity of audio is inversed analogously to experiment
for nidx=1:n_stim
    cal_audio = vertcat(cal_audio,vertcat(polarity*audio,zeros(jitter(nidx),1)));
    polarity  = -polarity;
end
clear audio

%% SoundmexPro
%--------------------------------------------------------------------------

% Apply Calibration
%------------------
% calibration for transient based on peak-value
level_stim = 20*log10(max(abs(cal_audio)));

% Compute initial gain to achieve start_val 
% level_stim + cal_val + gain_dB = start_val
% -> scaling = 10^(gain_dB/20)
gain_dB       = start_val - level_stim - cal_val;
current_level = start_val; % init
gain2add      = 0; % init

num_tracks = length(channel);
if num_tracks > 1
    tracks    = [0,1];
    plotcolor = {'b','r'};
else
    tracks = 0;
    if channel==0 % left
        gain_dB   = gain_dB(1);
        plotcolor = {'b'};
    end
    if channel==1 % right
        gain_dB   = gain_dB(2);
        plotcolor = {'r'};
    end
end

ok = soundmexpro('init', ...            % command name
    'force',      1, ...                % exit internally called before init
    'driver',     settings.asio_driver, ...
    'samplerate', settings.samp_freq, ...
    'output',     channel, ...
    'input',      [], ...
    'track',      num_tracks, ...
    'autocleardata', 1);            % clear audio data on memory
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

% Set variables to track state
%-----------------------------
max_runs       = 7;
direction      = -1; % go down
prev_direction = direction; 
correct_count  = 0; % number of correct answers

levels    = []; % store all levels
reversals = []; % store reversal levels for computation of hearing threshold
n_runs    = length(reversals); % A seriesof steps in one direction only is definedas a run
n_trials  = 0;

% Figure for plot
fig = figure;
movegui(fig,'northeast');

while n_runs < max_runs

    % Increment trial counter
    n_trials = n_trials +1;

    ok = soundmexpro('cleardata');
    if ~ok, error('error calling ''cleardata'' '); end

    pause(2)

    % Load signal
    %------------
    loopcount = 0; % play sounds once
    ok = soundmexpro('loadmem', ...     % command name
        'data', cal_audio, ...             % data vector
        'track',tracks, ...
        'loopcount', loopcount);
    if ~ok, error('error calling ''loadmem'''); end

    % Addjust gain
    %-------------
    % Signal Level in (p-p) peSPL
    current_level = current_level + gain2add;
    gain          = 10.^(gain_dB/20);

    % Signal levelin dB FS (peak)
    fprintf('Current level in dB FS (peak): %.1f\n', level_stim + gain_dB)

    % Check maximum level
    if current_level > settings.level_max
        error("Maximum level exceeded (%i > %i)!",current_level,settings.level_max)
        break
    end

    [ok, volume] = soundmexpro('trackvolume', ...    
        'track', tracks, ...
        'value', gain ...                           % volume values 
        );
    if ~ok, error('error calling ''trackvolume'''); end
    fprintf('Current level in dB FS (peak) via SoundMexPro: %.1f\n', 20*log10(volume))
    
    % Start SoundMexPro
    %------------------
    ok = soundmexpro('start','length',0); % device is never stopped, zeros are played endlessly 
    if ~ok,  error('error calling ''start''' ); end

    % Response - signal should be played until answer is given (or signal
    % ends - signal is played endlessly)
    %-----------------------------------
    title(sprintf('Current level: %.1f dB FS (peak) | Number of reversals = %i', current_level, n_runs),'FontSize',12);
    answer = questdlg('Do you hear the stimulus sequence?', ...
	         'Hearing threshold', ...
	         'Yes','No','Abort','Yes');

    switch answer
        case 'Yes'
            plotsymbol  = '+';
            response    = 1;
        case 'No'
            plotsymbol  = '_';
            response    = 0;
        case 'Abort'
            break % Ends execution of while loop -> soundmexpro is closed and results are stored
    end

    ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end 

    % Apply adaptive procedure and add gain
    %--------------------------------------
    [correct_count, direction, reversals, gain2add] = adaptive_staircase(response, correct_count, direction, reversals);

    gain_dB = gain_dB + gain2add; 

    % Check and append for reversals and count
    %-----------------------------------------
    if direction ~= prev_direction
        reversals = [reversals, current_level];
        n_runs    = length(reversals);
        disp(fprintf('Number of reversals: %i', n_runs))
    end

    % Plot Results
    %-------------
    % save data
    levels = horzcat(levels,vertcat(n_trials,response,current_level'));
    
    plot(levels(1,n_trials),levels(3,n_trials),plotsymbol,'MarkerSize',16,'LineWidth',4,'Color',plotcolor{1})
    hold on
    plot(levels(1,1:n_trials),levels(3,1:n_trials),'--','Color',plotcolor{1})
    if num_tracks>1
        plot(levels(1,n_trials),levels(4,n_trials),plotsymbol,'MarkerSize',16,'LineWidth',4,'Color',plotcolor{2})
        plot(levels(1,1:n_trials),levels(4,1:n_trials),'--','Color',plotcolor{2})
    end
    xlabel('Trial number','FontSize',10)
    ylabel('Level / dB (p-p) peSPL','FontSize',10)
end

ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end

disp('Playback of calibration_dBSL.m successfully completed.')

% Calculate the threshold as the mean of reversal levels
%-------------------------------------------------------
threshold_estimate = round(mean(reversals(3:end)),1);
title(sprintf('Threshold: %.1f dB (p-p) peSPL', threshold_estimate),'FontSize',16);
fprintf('Estimated Threshold: %.1f dB\n', threshold_estimate);

%% Save used Settings
%--------------------------------------------------------------------------

% Ask if you want to save the data
%---------------------------------

answer = questdlg("Should the measurement data be saved?",'Save','Yes','No','Yes');
switch answer
    case 'Yes'
        p   = fileparts(mfilename('fullpath'));
        pth = fullfile(p,'..','results',subjectname);
        if ~exist(pth,'dir')
            mkdir(pth)
        end
        results.Settings            = settings;
        results.Measdate            = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
        results.levels              = levels;
        results.reversals           = reversals;
        results.threshold_estimate  = threshold_estimate;
        results.used_channels       = channel;
        
        if isequal(channel,[0,1])
            channelname         = 'both';
            levels_legend = {'run';'heard (yes:1, no:0)';'level (left) / dB FS (RMS)';'level (right) / dB FS (RMS)'};
        elseif isequal(channel,0)
            channelname         = 'left';
            levels_legend = {'run';'heard (yes:1, no:0)';'level / dB FS (RMS)'};
        elseif isequal(channel,1)
            channelname         = 'right';
            levels_legend = {'run';'heard (yes:1, no:0)';'level / dB FS (RMS)'};
        end

        results.levels_legend = levels_legend;
        
        % Avoid overwriting of existing file
        in_file  = sprintf('%s_stim-%s_levels_%s_info.mat',subjectname,filename,channelname); 
        in_path  = pth;
        out_file = avoidOverwrite(in_file,in_path,2,0);
        save(fullfile(pth,out_file),'results')
    case 'No'   

end

close all

%% Extra functions
%--------------------------------------------------------------------------

    function [correct_count, direction, reversals, gain2add] = adaptive_staircase(response, correct_count, direction, reversals)
    %----------------------------------------------------------------------
    % ADAPTIVE_STAIRCASE Updates parameters in a 2-up, 1-down adaptive staircase procedure.
    %
    % Input:
    %   - response: Boolean (1 for correct, 0 for incorrect).
    %   - correct_count: Integer count of consecutive correct responses.
    %   - direction: Integer (1 for up, -1 for down).
    %   - reversals: Vector of levels at which reversals occurred.
    %
    % Output:
    %   - correct_count: Updated count of consecutive correct responses.
    %   - direction: Updated stimulus adjustment direction.
    %   - reversals: Updated vector of reversal levels.
    %   - gain2add: Gain adjustment for the stimulus level.
    %
    % Procedure:
    %   - Increment `correct_count` for correct responses.
    %   - After 3 correct responses, reset the count, set direction to down,
    %     and calculate `gain2add` using `gainrule`.
    %   - For incorrect responses, reset `correct_count`, set direction to up,
    %     and calculate `gain2add`.
    %
    %----------------------------------------------------------------------

        if response % correct answer
            correct_count = correct_count + 1;
            if correct_count == 2 % Go down after 3 correct responses
    
                % Reset counter
                correct_count = 0;
                % Remember direction for reversals
                % Go down
                direction = -1;
                % Adjust gain
                gain2add = gainrule(reversals, direction);
    
            else % Apply no gain
                gain2add = 0;
            end
    
        else % Incorrect answer
            correct_count = 0;  % Reset correct counter
            direction = 1;
            gain2add  = gainrule(reversals, direction);
        end

    end

    function gain2add = gainrule(reversals,direction)
    %----------------------------------------------------------------------
    % GAINRULE Determines the gain adjustment for the stimulus level.
    %
    % Input:
    %   - reversals: Vector of reversal levels.
    %   - direction: Integer indicating adjustment direction (-1 or 1).
    %
    % Output:
    %   - gain2add: Gain adjustment (positive or negative) based on the number
    %     of reversals and direction.
    %
    % Gain Adjustment Rules:
    %   - <2 reversals: 10 dB (large step)
    %   - 2-3 reversals: 5 dB (medium step)
    %   - ≥4 reversals: 2 dB (small step)
    %----------------------------------------------------------------------

        n_reversals = length(reversals);
        if n_reversals < 2
            gain2add = 10;
        elseif n_reversals == 2 || n_reversals == 3
            gain2add = 5;    
        elseif n_reversals >= 4
            gain2add = 2;
        end

        gain2add = direction*gain2add;

    end

end

