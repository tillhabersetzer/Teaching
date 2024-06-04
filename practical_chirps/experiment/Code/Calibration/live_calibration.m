function [rms_stim_dB] = live_calibration(channel, gaindB)
%--------------------------------------------------------------------------
% This function is used to measure the hearing threshold of a subject for
% a given signal. Here, the signal is a sequence of up-chirps with an
% ISI of 350ms + 50ms Jitter. This signal can be adjusted with a gain until
% it matches the hearing threshold. Positive values lead to an
% amplification, negative values to an attenuation. The signal can either
% be presented monaurally (left/right) or diotically. The level of the
% adjusted signal is returned. If the level matches the hearing threshold,
% than this level (rms_stim_dB) corresponds to 0 dB SL. Therefore, the
% CalVal is: CalVal = 0 dB - rms_stim_dB = -rms_stim_dB
% 
% output:
%   rms_stim_dB: Level in dB FS (RMS) of presented signal
%
% input:
%   channel: 0 left
%            1 right
%            [0,1] both
%   gaindB: channel gain in dB
%           >0: amplification
%           <0: attenuation
%--------------------------------------------------------------------------

% TO DO
% add limitation of gain: reasonable boundaries

%% Settings
%--------------------------------------------------------------------------

% Set default function paramters
if nargin<2
    gaindB = 0; % left
end
if nargin<1
    channel = 0; % left
end

% Import Calibration Signal
calibration_signal = 'calibration_up_rmseq.wav';
[calsig, SampFreq] = audioread(fullfile('..','Stimuli','calibration',calibration_signal));

% Import Settings
Settings   = importdata(fullfile('..','Settings','Settings.mat'));
ASIODriver = Settings.Driver;

% Other settings
loopcount = 0; % endless loop or 1 repitition

if Settings.SampFreq~=SampFreq
    error("Sampling frequencies don't match (%i~=%i).",Settings.SampFreq,SampFreq)
end

%% SoundmexPro
%--------------------------------------------------------------------------

% Attenuated calibration signal
% Pstim = AttdB + RMSdBStim
% -> scaling = 10^(AttdB/20)
Calsig = 10.^(gaindB/20)*calsig;

% Level of initial calibration signal in dB FS (RMS)
rms_stim_dB_init = 20*log10(rms(calsig));
% Level of scaled signal in dB FS (RMS)
rms_stim_dB = 20*log10(rms(Calsig));

disp(['Initial Stimulus Level: ', num2str(rms_stim_dB_init),' dB FS (RMS).'])
disp(['Stimulus Level is changed by ',num2str(gaindB), ' dB.'])
disp(['Adjusted stimulus Level: ',num2str(rms_stim_dB),' dB FS (RMS).'])
disp(['Note down Level of presented signal when hearing threshold is reached.'])

ok = soundmexpro('init', ...            % command name
    'force',      1, ...                % exit internally called before init
    'driver',     ASIODriver, ...
    'samplerate', SampFreq, ...
    'output',     channel, ...
    'input',      [], ...
    'autocleardata', 1);            % clear audio data on memory
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

ok = soundmexpro('loadmem', ...     % command name
    'data', Calsig, ...             % data vector
    'loopcount', loopcount);
if ~ok, error('error calling ''loadmem'''); end
   
ok = soundmexpro('start', ...
    'length',0); % device is never stopped, zeros are played endlessly
if ~ok,  error('error calling ''start''' ); end

answer = questdlg("Press 'Yes' to stop the playback.",'Input','Yes','Yes');
switch answer
    case 'Yes'
        ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
        ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end 

end
