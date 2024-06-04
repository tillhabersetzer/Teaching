function [rms_stim_dB] = calibration_dBSL(name,channel)
%--------------------------------------------------------------------------
% This function is used to measure the hearing threshold of a subject for
% a given signal adaptively. Here, the signal is a sequence of up-chirps 
% with an ISI of 350ms + 50ms Jitter. This signal is adjusted with a gain
% until the hearing threshold is found. The signal can either be presented 
% monaurally (left/right) or diotically. The level of the adjusted signal 
% is returned. If the level matches the hearing threshold, than this level 
% (rms_stim_dB) corresponds to 0 dB SL. Therefore, the CalVal is: 
% CalVal = 0 dB - rms_stim_dB = -rms_stim_dB
%
% The stepsize is 10 dB at first and is decreased to 5,3 and 1 dB when a
% different turning points have been passed.
% 
% output:
%   rms_stim_dB: Level in dB FS (RMS) of presented signal. Digital stimulus
%                level that corresponds to hearing threshold.
%   
%                Adaptively fitted curve with level values is saved in 
%                results folder.
%
% input:
%   name:    'subjectXX'
%   channel: 0 left
%            1 right
%            [0,1] both
%
%--------------------------------------------------------------------------

% TO DO
% add limitation of gain: reasonable boundaries
% check if trackvolume is applied appropriately
% set initial gain value -> set to dB SPL (peak)

%% Settings
%--------------------------------------------------------------------------
close all % close other figures

% Set default function paramters
if nargin<2
    channel = 0; % left
end
if nargin<1
    name = 'dummy';
end

% Import Calibration Signal
calibration_signal = 'calibration_up_rmseq.wav';
[calsig, SampFreq] = audioread(fullfile('..','Stimuli','calibration',calibration_signal));

% Import Settings
Settings   = importdata(fullfile('..','Settings','Settings.mat'));
ASIODriver = Settings.Driver;
% ASIODriver = 'Focusrite USB ASIO';
% ASIODriver = 'ASIO4ALL v2';

if Settings.SampFreq~=SampFreq
    error("Sampling frequencies don't match (%i~=%i).",Settings.SampFreq,SampFreq)
end

addpath(fullfile('..','Code','helper_functions'))

%% SoundmexPro
%--------------------------------------------------------------------------

% Attenuated calibration signal
% Pstim = AttdB + RMSdBStim
% -> scaling = 10^(AttdB/20)

% Start Value for gain (Find out!)
% GaindB = 0;
GaindB = Settings.Calibration.peakSPL.AttdB;

num_tracks = length(channel);
if num_tracks > 1
    tracks    = [0,1];
    plotcolor = {'b','r'};
else
    tracks = 0;
    if channel==0 % left
        GaindB    = GaindB(1);
        plotcolor = {'b'};
    end
    if channel==1 % right
        GaindB    = GaindB(2);
        plotcolor = {'r'};
    end
end

ok = soundmexpro('init', ...            % command name
    'force',      1, ...                % exit internally called before init
    'driver',     ASIODriver, ...
    'samplerate', SampFreq, ...
    'output',     channel, ...
    'input',      [], ...
    'track',      num_tracks, ...
    'autocleardata', 1);            % clear audio data on memory
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

% Counter for sign changes and runs
maxruns     = 20;
maxturns    = 8;
count_runs  = 0;
count_turns = 0;

% Store data
hearinglevel = []; % count_runs, response and rms_stim_dB

% Figure for plot
fig = figure;
movegui(fig,'northeast');

while count_turns < maxturns && count_runs<maxruns

    ok = soundmexpro('cleardata');
    if ~ok, error('error calling ''cleardata'' '); end

    pause(2)

    % Load signal
    %------------
    loopcount = 0; % play sounds once
    ok = soundmexpro('loadmem', ...     % command name
        'data', calsig, ...             % data vector
        'track',tracks, ...
        'loopcount', loopcount);
    if ~ok, error('error calling ''loadmem'''); end

    % Adjust gain
    %------------
    % Signal Level in dB FS (RMS)
    rms_stim_dB = 20*log10(rms(calsig)) + GaindB;
    gain        = 10.^(GaindB/20);

    [ok, ~] = soundmexpro('trackvolume', ...    
        'track', tracks, ...
        'value', gain ...                           % volume values 
        );
    if ~ok, error('error calling ''trackvolume'''); end
    
    % Start SoundMexPro
    %------------------
    ok = soundmexpro('start','length',0); % device is never stopped, zeros are played endlessly 
    if ~ok,  error('error calling ''start''' ); end

    % Response - signal should be played until answer is given (or signal
    % ends - signal is played endlessly)
    %-----------------------------------
    title(sprintf('Aktueller Pegel: %.2f dB FS (RMS)', rms_stim_dB(1)),'FontSize',16);
    answer = questdlg('Hören Sie die Stimulussequenz?', ...
	         'Hörschwelle', ...
	         'Ja','Nein','Abbruch','Ja');

    switch answer
        case 'Ja'
            gaindB = gainrule(count_turns,'Ja');
            plotsymbol = '+';
            response   = 1;
        case 'Nein'
            gaindB = gainrule(count_turns,'Nein');
            plotsymbol = '_';
            response   = 0;
        case 'Abbruch'
            break % ends execution of while loop -> soundmexpro is closed and results are stored
    end

    ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end   

    % Add Gains
    %----------
    GaindB = GaindB + gaindB; 

    % Set counters
    %-------------
    count_runs = count_runs+1;
    
    if count_runs > 1
        if sign(gaindB*gaindB_old)<0 % change in sign -> turning point
            count_turns = count_turns+1;
            disp(fprintf('Number of turns: %i', count_turns))
        end
    end
    gaindB_old = gaindB;

    % Plot Results
    %-------------
    % save data
    hearinglevel = horzcat(hearinglevel,vertcat(count_runs,response,rms_stim_dB'));
    

    plot(hearinglevel(1,count_runs),hearinglevel(3,count_runs),plotsymbol,'MarkerSize',16,'LineWidth',4,'Color',plotcolor{1})
    hold on
    plot(hearinglevel(1,1:count_runs),hearinglevel(3,1:count_runs),'--','Color',plotcolor{1})
    if num_tracks>1
        plot(hearinglevel(1,count_runs),hearinglevel(4,count_runs),plotsymbol,'MarkerSize',16,'LineWidth',4,'Color',plotcolor{2})
        plot(hearinglevel(1,1:count_runs),hearinglevel(4,1:count_runs),'--','Color',plotcolor{2})
    end
    xlabel('Darbietung','FontSize',10)
    ylabel('Pegel / dB FS (RSMS)','FontSize',10)
end

ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end

disp('Wiedergabe von calibration_dBSL.m erfolgreich beendet.')

%% Save used Settings
%--------------------------------------------------------------------------

% Ask if you want to save the data
%---------------------------------

answer = questdlg("Sollen die Messdaten gespeichert werden?",'Speichern','Ja','Nein','Ja');
switch answer
    case 'Ja'
        p   = fileparts(mfilename('fullpath'));
        pth = fullfile(p,'..','..','Results',name);
        if ~exist(pth,'dir')
            mkdir(pth)
        end
        info.Settings            = Settings;
        info.Measdate            = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
        info.hearinglevel        = hearinglevel;
        info.used_channels       = channel;
        
        if isequal(channel,[0,1])
            channelname         = 'both';
            hearinglevel_legend = {'run';'heard (yes:1, no:0)';'level (left) / dB FS (RMS)';'level (right) / dB FS (RMS)'};
        elseif isequal(channel,0)
            channelname         = 'left';
            hearinglevel_legend = {'run';'heard (yes:1, no:0)';'level / dB FS (RMS)'};
        elseif isequal(channel,1)
            channelname         = 'right';
            hearinglevel_legend = {'run';'heard (yes:1, no:0)';'level / dB FS (RMS)'};
        end

        info.hearinglevel_legend = hearinglevel_legend;
        
        % Avoid overwriting of existing file
        inFile  = [name,'_hearinglevel_',channelname,'_info.mat'];
        inPath  = pth;
        outFile = avoidOverwrite(inFile,inPath,2,0);
        save(fullfile(pth,outFile),'info')
    case 'Nein'     
end

%% extra functions
%--------------------------------------------------------------------------
    function gaindB = gainrule(count_turns,answer)
        if count_turns>4
            gaindB = 1; 
        elseif count_turns>3
            gaindB = 2; % 2dB
        elseif count_turns>2 
            gaindB = 5; % 5dB
        else
            gaindB = 10;
        end

        switch answer
            case 'Ja'
                gaindB = -gaindB; % attenuation
            case 'Nein'
                                  % amplification
        end
    end

end
