function experiment_click_basic()
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

% Load Stimulus
%==========================================================================
p               = fileparts(mfilename('fullpath'));
[stim,SampFreq] = audioread(fullfile(p,'..','..','Stimuli','rms_equalized','click_rmseq.wav'));
Settings        = importdata(fullfile(p,'..','..','Settings','Settings.mat'));

if Settings.SampFreq~=SampFreq
    error("Sampling frequencies don't match (%i~=%i).",Settings.SampFreq,SampFreq)
end

% Apply Calibration
%==========================================================================
% AttdB + RMSdBStim + CalVal = TargetLevel
% scaling = 10^(AttdB/20)

AttFact_stim = 10.^(Settings.Calibration.AttdB/20);

Stim      = zeros(length(stim),2);
Stim(:,1) = stim * AttFact_stim(1); % left
Stim(:,2) = stim * AttFact_stim(2); % right

% Generate Trigger
%==========================================================================

TrigID       = Settings.TrigID.click;                 
TrigLen      = round(SampFreq*0.01); % 10 ms
% TrigLen      = length(Stim);
TrigInfo     = EEGTrigID2info([TrigID],[16],'MEG'); % [8 8] distribution of 16 bits
TrigAmp      = TrigInfo.TrigWord;
stim_trigger = TrigAmp*ones(TrigLen,1);

Jitter = Settings.Jitter;
rng('shuffle')
Nstim          = Settings.Nstimuli;
% remove 10 ms because of extended stimulus
JitterDuration = (Jitter(1)-0.01) + Jitter(2)*rand(1,Nstim) ; % [0.34,0.39]s

% Extend Stimulus to 10 ms
%=========================
Stim = vertcat(Stim,zeros(length(stim_trigger)-length(Stim),2));

% Stimulus Presentation with SoundMexPro 
%==========================================================================
ASIODriver = Settings.Driver;

ok = soundmexpro('exit');
if ~ok, error('error calling ''exit'' '); end             

ok = soundmexpro('init', ...       % command name
    'samplerate',   SampFreq, ...
    'force',        1, ...         % exit called internally before init
    'driver',       ASIODriver, ...
    'output',       [0 1 8], ...   % soundcard channels [left,right,trigger]
    'input',        -1, ...        % no input channel used
    'track',         3, ...        % [stim stim trigger] [left right trigger]
    'autocleardata', 1);           % audio data already played completely should be cleared from memory automatically on next data loading command
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('trackmap','track', [0 1 2]); % virtual channel (first track -> channel 0, second track -> channel 1, third track -> channel 3)        
if ~ok, error('error calling ''trackmap'' '); end             
 
ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

% Play Zeros
%==========================================================================

% switch on debug saving
if 1 ~= soundmexpro('debugsave', ...    % command name
        'value', [1], ...               % enable/disable flag (here: enable)
        'channel', [2] ...              % channel(s) where to apply value
        );
    error(['error calling ''debugsave''' error_loc(dbstack)]);
end

ok = soundmexpro('loadmem', ...
     'data',zeros(1*SampFreq,3), ...
     'track',[0 1 2], ...
     'loopcount', 1);
if ~ok,  error('error calling ''loadmem''' ); end

ok = soundmexpro('start', ...
    'length', 0);              % device is never stopped, zeros are played endlessly. In this case you may load new data to track(s) at any time.
if ~ok,  error('error calling ''start''' ); end

pause(1)
    
% Play Stimuli
%==========================================================================

% progress bar
fig = uifigure;
d = uiprogressdlg(fig,...
                  'Title','Fortschritt',...
                  'Icon','info',...
                  'ShowPercentage','on',...
                  'Cancelable','on');

% Deliver 50% of stimulus with reversed phase
flipsign = zeros(1,Nstim);
fsign    = 1;

for n = 1:Nstim
    
    if d.CancelRequested
             break
    end
    
    ok = soundmexpro('loadmem', ...            % command name
        'data', [fsign*Stim, stim_trigger], ...      % matrix with one or more columns of data
        'track',[0 1 2], ...                   % vector with tracks, where data to be played
        'loopcount', 1);                       % 0 is an endless loop
    if ~ok,  error('error calling ''loadmem'' '); end
    
    pause(length(Stim)/SampFreq + JitterDuration(n)); 
    d.Value   = n/Nstim;
    d.Message =['Stimulus ',num2str(n),' von ',num2str(Nstim),' gespielt.'];  

    if d.CancelRequested
         break
    end
    
    % flip stimulus
    flipsign(n) = fsign;
    fsign       = -fsign;
end

ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end   

disp('Wiedergabe von experiment_click_basic.m erfolgreich beendet.')

end
