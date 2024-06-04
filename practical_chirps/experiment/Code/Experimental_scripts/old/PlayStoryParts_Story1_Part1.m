function PlayStoryParts_Story1_Part1(f,Settings)
%==========================================================================
% plays story part
%
% parameter: f figure handle for progress bar
%            Settings: file with necessary information 
%            - last update
%            - driver
%            - calibration values
%==========================================================================

% Load Files
%==========================================================================

% check_track = 0; % option to check loaded data on tracks with soundmexpro

p = fileparts(mfilename('fullpath'));
[story, SampFreq] = audioread([p filesep '..' filesep 'Geschichte' filesep 'Das_Schwatzende_Herz' filesep 'Das_Schwatzende_Herz1_final.wav']);  
[tone, SampFreq2] = audioread([p filesep '..' filesep 'ton_stimuli' filesep Settings.tone '.wav']);
[start,SampFreq3] = audioread([p filesep '..' filesep 'ton_stimuli' filesep 'Start.wav']);
[ende,SampFreq4]  = audioread([p filesep '..' filesep 'ton_stimuli' filesep 'End.wav']);

if ~isequal(SampFreq,SampFreq2,SampFreq3,SampFreq4)
    error(['SampFreq of loaded data are not the same'])
end

% Apply Calibration
%==========================================================================
% story CalVal
CalVal     = Settings.CalVal;
StoryLevel = Settings.StoryLevel;

RMSdBstory = 20.*log10(rms(story));
AttdB      = (CalVal-StoryLevel+RMSdBstory);
AttFact    = 10.^(-AttdB/20);  

% tone CalVal
CalVal_tone = Settings.CalVal_tone;
ToneLevel   = Settings.ToneLevel;

RMSdBtone    = 20.*log10(rms(tone));
AttdB_tone   = (CalVal_tone-ToneLevel+RMSdBtone);
AttFact_tone = 10.^(-AttdB_tone/20);  

Story      = zeros(length(story),2);
Story(:,1) = story * AttFact(1);
Story(:,2) = story * AttFact(2);

Tone      = zeros(length(tone),2);
Tone(:,1) = tone * AttFact_tone(1);
Tone(:,2) = tone * AttFact_tone(2);

Start      = zeros(length(start),2);
Start(:,1) = start * AttFact(1);
Start(:,2) = start * AttFact(2);

End      = zeros(length(ende),2);
End(:,1) = ende * AttFact(1);
End(:,2) = ende * AttFact(2);

% Generate Trigger
%==========================================================================
% story
TrigWin                  = SampFreq*1;                     % trigger every 1 seconds
TrigLen                  = 0.2*SampFreq;
TrigInfo                 = EEGTrigID2info([1],[16],'MEG'); % [8 8] distribution of 16 bits
TrigAmp                  = TrigInfo.TrigWord;    
story_trigger            = zeros(TrigWin,1);
story_trigger(1:TrigLen) = TrigAmp;
L                        = length(Story);

% tone                   
TrigLen      = length(tone);                   % trigger every tone pulse             
TrigInfo     = EEGTrigID2info([2],[16],'MEG'); % [8 8] distribution of 16 bits
TrigAmp      = TrigInfo.TrigWord;
tone_trigger = TrigAmp*ones(TrigLen,1);

rng('shuffle')
Ntones         = Settings.Ntones;
JitterDuration = 1.1 + 0.2*rand(1,Ntones) ; % [1.1,1.3]

% Get SoundMexPro Into Work
%==========================================================================
ok = soundmexpro('exit');
if ~ok, error('error calling ''exit'' '); end             

ASIODriver = Settings.Driver;

ok = soundmexpro('init', ...     % command name
    'samplerate', SampFreq, ...
    'force',      1, ...         % exit called internally before init
    'driver',     ASIODriver, ...
    'output',     [0 1 8], ...   % soundcard channels
    'input',      -1, ...        % no input channel used
    'track',       3);           % [story story trigger] [links rechts trigger]
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('trackmap','track', [0 1 2]); % virtual channel          
if ~ok, error('error calling ''trackmap'' '); end             
 
ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

%==========================================================================
% if check_track
% if 1 ~= soundmexpro('showtracks', 'wavedata', 0, 'topmost', 0)
%     error(['error calling ''showtracks''' error_loc(dbstack)]);
% end
% end

% Play Zeros
%==========================================================================
ok = soundmexpro('start','length',0);
if ~ok,  error(['error calling ''start''' ]); end

% Load Start
%==========================================================================
ok = soundmexpro('loadmem', ...
     'data',Start, ...
     'track',[0 1], ...
     'loopcount', 1);
if ~ok,  error(['error calling ''loadmem''' ]); end
 
% if check_track
% ok = soundmexpro('updatetracks');
% if ~ok,  error(['error calling ''updatetracks'' ']); end
% end

% wait
pause(length(Start)/SampFreq + 2);
    
% Play Tones
%==========================================================================
nbytes = fprintf('tone 0 of %d', Ntones);
for n = 1:Ntones
    
    fprintf(repmat('\b',1,nbytes))
    nbytes = fprintf('tone %d of %d\n', n, Ntones);
    
    ok = soundmexpro('loadmem', ...            % command name
        'data', [Tone, tone_trigger], ...      % data vector
        'track',[0 1 2], ...
        'loopcount', 1);
    if ~ok,  error(['error calling ''loadmem'' ']); end
    
%     if check_track
%     ok = soundmexpro('updatetracks');
%     if ~ok,  error(['error calling ''updatetracks'' ']); end
%     end
    
    pause(length(tone)/SampFreq + JitterDuration(n));  
    
end

ok = soundmexpro('stop');
if ~ok, error('error calling ''stop'' '); end

% wait 4 sec
pause(4)

% Play Story
%==========================================================================

% add 1 sec zeros first to avoid wrongly setted triggers
ok = soundmexpro('loadmem', ...        % command name
    'data', zeros(SampFreq,3), ...       % data vector
    'track',[0,1,2], ...
    'loopcount', 1);
if ~ok,  error(['error calling ''loadmem'' ']); end

ok = soundmexpro('loadmem', ...        % command name
    'data', [story_trigger], ...       % data vector
    'track',[2], ...
    'loopcount', 0);
if ~ok,  error(['error calling ''loadmem'' ']); end

ok = soundmexpro('loadmem', ...      % command name
    'data', [Story], ...             % data vector
    'track',[0 1], ...
    'loopcount', 1);
if ~ok,  error(['error calling ''loadmem'' ']); end

% if check_track
% ok = soundmexpro('updatetracks');
% if ~ok,  error(['error calling ''updatetracks'' ']); end
% end

ok = soundmexpro('start','length', L); % length, here -1: stop if no track has more data (this is default));         
if ~ok,  error(['error calling ''start'' ']); end

basic_progress_bar(f,L,SampFreq,5) % update every 5s

soundmexpro('stop');

% Load End
%==========================================================================
% wait 2 sec
pause(2);

ok = soundmexpro('loadmem', ...
     'data',End, ...
     'track',[0 1], ...
     'loopcount', 1);
if ~ok,  error(['error calling ''loadmem''' ]); end

ok = soundmexpro('start','length',0);
if ~ok,  error(['error calling ''start''' ]); end

ok = soundmexpro('wait','track',0,'mode','output');
if ~ok,  error(['error calling ''wait'' ']); end

soundmexpro('stop');
soundmexpro('exit');

disp("Wiedergabe von PlayStoryParts_Story1_Part1.m erfolgreich beendet")

end