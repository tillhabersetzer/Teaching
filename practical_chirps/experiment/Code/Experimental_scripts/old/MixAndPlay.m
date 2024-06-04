function MixAndPlay(f,SRT,name,Settings,Sentences_data)
%--------------------------------------------------------------------------
% parameter: 

% f: figure handle for progress bar
% SRT: speech recognition threshold 50%
% name: name of participant

% Settings: file with necessary information 
% - last update
% - driver
% - calibration values

% Sentences_data
% - Sentences_data.Sentences         = sentences;
% - Sentences_data.TestList          = TestList;
% - Sentences_data.Intelligibility   = intelligibility;
% - Sentences_data.TriggerAssignment = TriggerAssignment; contains Settings
% about Trigger Assignments
%--------------------------------------------------------------------------

% check_track = 0; % option to check loaded data on tracks with soundmexpro

p = fileparts(mfilename('fullpath'));

CalVal     = Settings.CalVal; % dB bei virtueller Vollaussteurung von 0 dB FS (linker/rechter Kopfhörer)
NoiseLevel = 65; % Rauschpegel in dB

TestList        = Sentences_data.TestList;
intelligibility = Sentences_data.Intelligibility;

TestListSN = get_snr(SRT,intelligibility);

[noise,SampFreq]  = audioread([p filesep '..' filesep 'Noise' filesep 'olnoise_long.wav']);
[start,SampFreq2] = audioread([p filesep '..' filesep 'ton_stimuli' filesep 'Start.wav']);
[ende,SampFreq3]  = audioread([p filesep '..' filesep 'ton_stimuli' filesep 'End.wav']);

if ~isequal(SampFreq,SampFreq2,SampFreq3)
    error(['SampFreq of loaded data are not the same'])
end

sentences_data    = Sentences_data.Sentences;
TriggerAssignment = Sentences_data.TriggerAssignment;

rng('shuffle'); 
% die ersten beiden Sätze sollten verständlich sein, sonst weiß der Proband
% nicht wann die Messung beginnt (feedback Jasper)
% Liste 30 mit 95%-Verständlichkeit kommt zuerst, daher flipud, dann
% permutiere die Sätze der ersten Liste, da die ersten beiden Elemente der
% Playmatrix willkürlich aus Liste 30 kommen sollen. Dann permutiere noch
% einmal die verbleibenden 118 Sätze auf den Positionen 3-120
TriggerAssignment         = flipud(TriggerAssignment);
order1                    = randperm(20)';
TriggerAssignment(1:20,:) = TriggerAssignment(order1,:);

order2                     = (randperm(118)+2)';
TriggerAssignment(3:120,:) = TriggerAssignment(order2,:);

% TrigerAssignment liefert neue Playlist ->  mit Zeilen (TriggerID,Liste,Satz)

% Apply Calibration
%==========================================================================

TotalLenNoise = length(noise);
rng('shuffle')
act_noise   = circshift(noise,randi(TotalLenNoise));
act_noise   = act_noise(1:7*60*SampFreq); % only take the first 7 min
RMSdB_Noise = 20*log10(rms(noise));
AttdB       = (CalVal-NoiseLevel+RMSdB_Noise);
AttFact     = 10.^(-AttdB/20);  

Act_noise      = zeros(length(act_noise),2);
Act_noise(:,1) = act_noise * AttFact(1); % linker Kopfhörer
Act_noise(:,2) = act_noise * AttFact(2); % rechter Kopfhörer

Start      = zeros(length(start),2);
Start(:,1) = start * AttFact(1);
Start(:,2) = start * AttFact(2);

Ende      = zeros(length(ende),2);
Ende(:,1) = ende * AttFact(1);
Ende(:,2) = ende * AttFact(2);

% Get SoundMexPro Into Work
%==========================================================================
ok = soundmexpro('exit');
if ~ok, error('error calling ''exit'' '); end   

ASIODriver = Settings.Driver; 
TrigLen    = 0.2.*SampFreq;

ok = soundmexpro('init', ...     % command name
    'samplerate', SampFreq, ...
    'force',      1, ...         % exit called internally before init
    'driver',     ASIODriver, ...
    'output',     [0 1 8], ...
    'input',        -1, ...      % no input channel used
    'track',       5);           % [noise,signal,noise,signal,trigger] [links links recht rechts trigger]
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('trackmap','track', [0 0 1 1 2]);
if ~ok, error('error calling ''trackmap'' '); end          

ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

%==========================================================================
% if check_track
% if 1 ~= soundmexpro('showtracks', 'wavedata', 0, 'topmost', 0)
%     error(['error calling ''showtracks''' error_loc(dbstack)]);
% end
% end
% Load Start
%==========================================================================
ok = soundmexpro('loadmem', ...
     'data',[Start;zeros(2*SampFreq,2)], ...
     'track',[0 2], ...
     'loopcount', 1);
if ~ok,  error(['error calling ''loadmem''' ]); end

% Load Noise
%==========================================================================
ok = soundmexpro('loadmem', ...  % command name
    'data', Act_noise, ...       % data vector
    'track',[0 2], ...
    'loopcount', 1);
if ~ok,  error(['error calling ''loadmem''' ]); end

% if check_track
% ok = soundmexpro('updatetracks','wavedata',0);
% if ~ok,  error(['error calling ''updatetracks'' ']); end
% end

% Play Zeros
%==========================================================================
ok = soundmexpro('start','length',0);
if ~ok,  error(['error calling ''start''' ]); end

% wait 3 sec - avoid onset effects of noise in first olsa recording
pause(length([Start;zeros(2*SampFreq,2)])/SampFreq+3)

% Start Soundfiles
%==========================================================================

% progress bar
d = uiprogressdlg(f,...
                  'Title','Fortschritt',...
                  'Icon','info',...
                  'ShowPercentage','on',...
                  'Cancelable','on');

for k = 1:length(TriggerAssignment)
    
    if d.CancelRequested
             break
    end
    
    Lidx   = find(TriggerAssignment(k,2) == TestList); % Gibt Index im Vektor Testlist aus, also 1,2,3,4,5,6
    Sidx   = TriggerAssignment(k,3);    
    act_SN = TestListSN(Lidx); % feste Zuordnung zwischen SNR und Testliste
    
    JitterDuration    = 1 + 0.2.*rand(1,1); % whatever you want [1:1.2]
    act_sentence      = sentences_data{Lidx,Sidx};
    
    SentLevel         = NoiseLevel + act_SN;
    RMSdBSent         = 20.*log10(rms(act_sentence(:,1)));
    AttdB             = (CalVal-SentLevel+RMSdBSent);
    AttFact           = 10.^(-AttdB/20);  
    
    Act_sentence(:,1) = act_sentence(:,1) * AttFact(1); % sentence ist Stereo Signal (sind aber leicht unterschiedlich)
    Act_sentence(:,2) = act_sentence(:,1) * AttFact(2); % -> verwende nur einen Kanal davon

    SentLen           = length(act_sentence);
    
    TrigInfo          = EEGTrigID2info([TriggerAssignment(k,1)],[16],'MEG'); % [8 8] distribution of 16 bits
    TrigAmp           = TrigInfo.TrigWord;
    
    act_trigger            = zeros(SentLen,1);
    act_trigger(1:TrigLen) = TrigAmp;
    
    ok = soundmexpro('loadmem', ...                       % command name
        'data', [Act_sentence, act_trigger], ...          % data vector
        'track',[1 3 4], ...
        'loopcount', 1);
    
%     if check_track
%     ok = soundmexpro('updatetracks','wavedata',0);
%     if ~ok,  error(['error calling ''updatetracks'' ']); end
%     end
    
    pause(SentLen/SampFreq + JitterDuration);
    
    d.Value   = k/length(TriggerAssignment);
    d.Message =['Satz ',num2str(k),' von ',num2str(length(TriggerAssignment)),' gespielt.'];
    
    Act_sentence = [];
    
end

% wait 3 sec afterwards, so that there is sufficiently long noise after the last sentence 
pause(3);

soundmexpro('stop');

% wait 2 sec after noise offset
pause(2);

% Load End
%==========================================================================
ok = soundmexpro('loadmem', ...
     'data',Ende, ...
     'track',[0 2], ...
     'loopcount', 1);
if ~ok,  error(['error calling ''loadmem''' ]); end

ok = soundmexpro('start','length',0);
if ~ok,  error(['error calling ''start''' ]); end

ok = soundmexpro('wait','track',0,'mode','output');
if ~ok,  error(['error calling ''wait'' ']); end

soundmexpro('stop');
soundmexpro('exit');

% save parameters and informations in info struct
%==========================================================================

Table_PlayList = array2table(TriggerAssignment,...
    'VariableNames',{'TriggerID','Testliste','Satznummer'});
Table_Parameter = array2table([TestList',intelligibility',TestListSN'],...
    'VariableNames',{'Testliste','Verstaendlichkeit','SNR'});

info.Name            = name;
info.Date            = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')); % letztes Update der Einstellungen
info.SRT             = SRT;

info.Table_PlayList  = Table_PlayList;
info.Table_Paramater = Table_Parameter;
save([p filesep '..' filesep 'Results', filesep, name,'_info.mat'],'info')

disp("Wiedergabe von MixAndPlay.m erfolgreich beendet")

end
