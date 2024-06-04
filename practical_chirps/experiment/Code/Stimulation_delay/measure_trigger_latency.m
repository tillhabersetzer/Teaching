function measure_trigger_latency(channel,gaindB,stimtype)
%--------------------------------------------------------------------------
%
% input:
%   channel: 0 left (default)
%            1 right
%            
%   gaindB: channel gain in dB (default: 0)
%           <0: attenuation
%           >0: amplification (does not make sense, signals are already full scale)
%
%   stimtype: 'click' (default)
%             'up'
%             'down'
%--------------------------------------------------------------------------

% add limitation of gain: reasonable boundaries

%% Settings
%==========================================================================

% Set default function paramters
if nargin <3
%     stimtype = 'click';
    stimtype = 'click';
end
if nargin<2
    gaindB = 0; % left
end
if nargin<1
    channel = 0; % left
end

disp('hey')

% Triggerbox functions
addpath(fullfile('..','EEGTriggerTB'))

% Import Settings
Settings        = importdata(fullfile('..','Settings','Settings.mat'));
ASIODriver      = Settings.Driver;
SampFreq        = Settings.SampFreq;
Nstim           = 1000;
trigger_channel = 8;
% Focusrite
% ASIODriver = 'Focusrite USB ASIO';
% trigger_channel = 1;

% Jitter
jitterstat = 1; % 0 or 1

%% Signal Generation
%==========================================================================

% Repeated Stimulus Pattern
% [calsig, TrigLen-length(calsig)] + [Pause]
% --------------------------------
%       length of trigger
%
% -> pause is: TriLen - length(calsig) + Pause
% e.g.   20 ms - 10 ms + 350 ms = 360 ms

% Trigger 
TrigID   = 1; 
TrigLen  = round(SampFreq*0.02); % 20 ms >= length(calsignal)!
TrigInfo = EEGTrigID2info(TrigID,16,'MEG');
TrigAmp  = TrigInfo.TrigWord;
trigger  = TrigAmp*ones(TrigLen,1);

% Pause 
if jitterstat
    JitterDuration_sec = 0.35 + 0.05*rand(1,Nstim); % [0.35,0.40]s (Jitter)
    JitterDuration     = round(JitterDuration_sec*SampFreq); % samples
else
    JitterDuration     = round(0.375*SampFreq)*ones(1,Nstim); % samples 
end

switch stimtype
    case 'click'
        % ~ 4 samples duration with 44.1 kHz
        [calsig, SampFreq2] = audioread(fullfile('..','Stimuli','raw','click.wav'));
    case 'up'
        % ~ 10 ms duration
        [calsig, SampFreq2] = audioread(fullfile('..','Stimuli','raw','up.wav'));
    case 'down'
        % ~ 10 ms duration
        [calsig, SampFreq2] = audioread(fullfile('..','Stimuli','raw','down.wav'));
end

% Extend Stimulus to TrigLen for Alignment
calsig = vertcat(calsig,zeros(length(trigger)-length(calsig),1));

if SampFreq~=SampFreq
    error("Sampling frequencies don't match (%i~=%i).",Settings.SampFreq,SampFreq)
end

%% Signal Scaling and SoundMexPro
%==========================================================================

% Attenuated calibration signal
% Pstim = AttdB + RMSdBStim
% -> scaling = 10^(AttdB/20)
Calsig = 10.^(gaindB/20)*calsig;

% Level of initial calibration signal in dB FS (Peak)
rms_stim_dB_init = 20*log10(max(abs(calsig)));
% Level of scaled signal in dB FS (Peak)
rms_stim_dB = 20*log10(max(abs(Calsig)));

disp(['Initial Stimulus Level: ', num2str(rms_stim_dB_init),' dB FS (Peak).'])
disp(['Stimulus Level is changed by ',num2str(gaindB), ' dB.'])
disp(['Adjusted stimulus Level: ',num2str(rms_stim_dB),' dB FS (Peak).'])

ok = soundmexpro('init', ...            % command name
    'force',        1, ...              % exit internally called before init
    'driver',       ASIODriver, ...
    'samplerate',   SampFreq, ...
    'output',       [channel trigger_channel], ...
    'input',        -1, ...             % no input channels
    'track',         2,...
    'autocleardata', 1);                % clear audio data on memory
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('trackmap','track', [0 1]);     
if ~ok, error('error calling ''trackmap'' '); end 

ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

%% Play Stimuli
%==========================================================================

% Counter
tracks_loaded = 0;

% progress bar
fig = uifigure;
d = uiprogressdlg(fig,...
                  'Title','Signals loaded',...
                  'Icon','info',...
                  'ShowPercentage','on',...
                  'Cancelable','on');

ok = soundmexpro('start','length', 0);             
if ~ok,  error('error calling ''start''' ); end

while tracks_loaded<Nstim

    [ok, trackload] = soundmexpro('trackload');
    if ok ~= 1
        error('error calling ''trackload''');
        break;
    end

    if d.CancelRequested
         break
    end
    
    if (trackload(1) < 10)
        
       ok = soundmexpro('loadmem', ...    
            'data', vertcat([Calsig, trigger],zeros(JitterDuration(tracks_loaded+1),2)), ...             
            'track', [0 1], ...
            'loopcount', 1); ...
        if ~ok, error('error calling ''loadmem'''); end                    
                if ~ok  
                    error('error calling ''loadmem'' ');
                    break;
                end
                tracks_loaded  = tracks_loaded + 1;
                d.Value        = tracks_loaded/Nstim;
                d.Message      = ['Signal ',num2str(tracks_loaded),' of ',num2str(Nstim),' loaded.'];  
                pause(0.01)
    end 

end

pause(5) % pause is needed for preloaded buffer: 10 files
disp('Playback finished.')

ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end 

end
