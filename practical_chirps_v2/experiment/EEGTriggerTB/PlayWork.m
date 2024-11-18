

 %[TrigInfo] = EEGTrigID2info( {[1 2 3]},'idx',24,8);
 [TrigInfo] = EEGTrigID2info( '00001111','bitmask',24,8)
 [TrigInfo] = EEGTrigID2info( '10000000','bitmask',24,8)
 [TrigInfo] = EEGTrigID2info( '11110000','bitmask',24,8)
 [TrigInfo] = EEGTrigID2info( '00000010','bitmask',24,8)
 
 
  TrigInfo.Trigword

 
 
ASIODriver = 'ASIO Fireface USB';
SampFreq   = 44100;

%---------------------------------------
ok = soundmexpro('exit');
%----
ok = soundmexpro('init', ...     % command name
    'force',      1, ...
    'driver',     ASIODriver, ...
    'output',     [0 1 8]);

if ~ok, error('MATLAB:MRTPLAYREC:CannotInitializeSoundmex','error calling ''init'''); end
%----
ok = soundmexpro('cleardata');
if ~ok, error('MATLAB:MRTPLAYREC:CannotCleardata','error calling ''cleardata'''); end
%----

ok = soundmexpro('start','length',0);
if ~ok,  error(['error calling ''start''' ]); end

%===========================================
 PlayData = zeros(500,3);
 PlayData(10:510,3) = TrigInfo.Trigword; 

 
ok = soundmexpro('loadmem', ...      % command name
    'data', PlayData, ...          % data vector
    'loopcount', 1);
if  ~ok, error('error calling ''loadmem'''); end
%===========================================


ok = soundmexpro('stop');
pause(1)
soundmexpro('exit')
