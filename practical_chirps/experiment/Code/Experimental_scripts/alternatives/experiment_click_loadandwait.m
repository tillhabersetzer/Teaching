function experiment_click_loadandwait()
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
JitterDuration = round(JitterDuration*SampFreq); % samples

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

% Play Stimuli
%==========================================================================

% progress bar
fig = uifigure;
d = uiprogressdlg(fig,...
                  'Title','Tracks loaded',...
                  'Icon','info',...
                  'ShowPercentage','on',...
                  'Cancelable','on');

% Length
L = sum(JitterDuration) + Nstim*length(Stim) + SampFreq*2;

% Counter
tracks_loaded = 0;

% Preload data
fsign = 1;
for n = 1:4
    ok = soundmexpro('loadmem', ...
     'data',vertcat([fsign*Stim, stim_trigger],zeros(JitterDuration(tracks_loaded+1),3)), ...
     'track',[0 1 2], ...
     'loopcount', 1);
    if ~ok,  error('error calling ''loadmem''' ); end

    tracks_loaded = tracks_loaded + 1;
    d.Value       = tracks_loaded/Nstim;
    d.Message     = ['Track ',num2str(tracks_loaded),' von ',num2str(Nstim),' geladen.']; 
    fsign         = -fsign;
    pause(0.01)
end

ok = soundmexpro('start','length', L);             
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
    
    if (trackload(1) < 4)
        
        ok = soundmexpro('loadmem', ...            
            'data', vertcat([fsign*Stim, stim_trigger],zeros(JitterDuration(tracks_loaded+1),3)), ...     
            'track',[0 1 2], ...                   
            'loopcount', 1);                       
        if ~ok  
            error('error calling ''loadmem'' ');
            break;
        end
        tracks_loaded = tracks_loaded + 1;
        d.Value       = tracks_loaded/Nstim;
        d.Message     = ['Track ',num2str(tracks_loaded),' von ',num2str(Nstim),' geladen.'];  
        fsign         = -fsign;
        pause(0.01)
    end 

end

% ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
% ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end   

disp('Wiedergabe von experiment_click_loadandwait.m erfolgreich beendet.')

end