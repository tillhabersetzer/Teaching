close all
clear all
clc

addpath(fullfile('..','..','code','EEGTriggerTB'))

TrigID       = 2;               
TrigLen      = round(44100*0.1);
TrigInfo     = EEGTrigID2info([TrigID],[16],'MEG'); % [8 8] distribution of 16 bits
TrigAmp      = TrigInfo.TrigWord;
stim_trigger = TrigAmp*ones(TrigLen,1);

stim_trigger = vertcat(stim_trigger,zeros(10000,1));


ok = soundmexpro('init', ...            % command name
    'force',      1, ...                % exit internally called before init
    'driver',     'ASIO Fireface USB', ...
    'samplerate', 44100, ...
    'output',     8, ...
    'input',      [], ...
    'autocleardata', 1);            % clear audio data on memory
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

ok = soundmexpro('loadmem', ...     % command name
    'data', stim_trigger, ...             % data vector
    'loopcount', 0);
if ~ok, error('error calling ''loadmem'''); end
   
ok = soundmexpro('start', ...
    'length',0); % device is never stopped, zeros are played endlessly
if ~ok,  error('error calling ''start''' ); end


% or use test mode of EEG Trigger-Toolbox
%--------------------------------------------------------------------------

addpath(fullfile('..','..','code','EEGTriggerTB'))

TrigID         = 2;   
BitSplitScheme = 16;
TrigDur        = 100; % ms

ASIO.Driver    = 'ASIO Fireface USB';
ASIO.SampFreq  = 44100;
ASIO.TrigChan  = 8;

TrigInfo       = EEGTrigID2info(TrigID,BitSplitScheme,'MEG',TrigDur,ASIO); % playse trigger with 10 ms pause in between

