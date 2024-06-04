function experiment_down(f,Settings,name)
%--------------------------------------------------------------------------
% parameters: 
%
% f:        figure handle
% Settings: experimental settings 
% name:     subject code
%--------------------------------------------------------------------------

% check_tracks = 0; % option to check loaded data on tracks with soundmexpro

% Load Stimulus
%==========================================================================
p               = fileparts(mfilename('fullpath'));
[stim,SampFreq] = audioread(fullfile(p,'..','Stimuli','rms_equalized','down_rmseq.wav'));

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

TrigID       = Settings.TrigID.downchirp;                 
%TrigLen      = round(SampFreq*0.01); % 10 ms
TrigLen      = length(Stim);
TrigInfo     = EEGTrigID2info([TrigID],[16],'MEG'); % [8 8] distribution of 16 bits
TrigAmp      = TrigInfo.TrigWord;
stim_trigger = TrigAmp*ones(TrigLen,1);

Jitter = Settings.Jitter;
rng('shuffle')
Nstim          = Settings.Nstimuli;
JitterDuration = Jitter(1) + Jitter(2)*rand(1,Nstim) ; % [0.35,0.4]s

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

% Check tracks
%==========================================================================
% shows visualization of files/vectors in tracks. This 'view' is
% especially intended to check the setup of your experiment/pardigm,
% i.e. if everything is loaded/located as expected. During playback
% a cursor shows the current position.

% if check_tracks
%     if 1 ~= soundmexpro('showtracks', ...
%             'wavedata', 0, ...            % if set to '1' waveforms are painted as well
%             'topmost', 0)                 % if set to '1' track window stays on top
%         error('error calling ''showtracks'''); 
%     end
% end

% Play Zeros
%==========================================================================

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
d = uiprogressdlg(f,...
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

%     if check_tracks     
%         ok = soundmexpro('updatetracks');
%         if ~ok,  error('error calling ''updatetracks'' '); end    
%     end
    
    pause(length(Stim)/SampFreq + JitterDuration(n)); 
    d.Value   = n/Nstim;
    d.Message =['Stimulus ',num2str(n),' von ',num2str(Nstim),' gespielt.'];  
    
    % flip stimulus
    flipsign(n) = fsign;
    fsign       = -fsign;
end

ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end   

disp('Wiedergabe von experiment_down.m erfolgreich beendet.')

% Save used Settings
%--------------------------------------------------------------------------
answer = questdlg("Sollen die Messdaten gespeichert werden?",'Speichern','Ja','Nein','Ja');
switch answer
    case 'Ja'
        pth = fullfile(p,'..','..','Results',name);
        if ~exist(pth,'dir')
            mkdir(pth)
        end
        info.Settings = Settings;
        info.Measdate = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')); % letztes Update der Einstellungen
        info.Jitter   = JitterDuration;
        info.flipsign = flipsign;
        
        % Avoid overwriting of existing file
        inFile  = [name,'_down_info.mat'];
        inPath  = pth;
        outFile = avoidOverwrite(inFile,inPath,2,0);
        save(fullfile(pth,outFile),'info')
    case 'Nein'
end

end

% Notes
%--------------------------------------------------------------------------

% loadfile vs loadmem
%---------------------
% Loadfile does not load the complete file to
% memory, therefore it should be used rather than 'loadmem'
% on huge files. But if files are small and to be played in
% loop it is recommended to use 'wavread' and 'loadmem' to
% play data, because it is much more efficient to read from
% memory than to read from file on the fly.