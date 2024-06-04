function [status,trigmode,TrigIDPlayed] = EEGTrigPlay(Mode, EEGPlayMatrix, Trig,  NumOfLoops, BlockOrder)

persistent status
persistent Trig
JitterInBlockSamp
JitterOutBlockSamp
LatencySamp

% Trig.ID
% Trig.IDMode
% Trig.BitRef
% Trig.NumOfBits,
% Trig.TwoBitScheme
% Trig.DurSamp
% Trig.LatSamp
% 
if isempty(status), status = [0 0 0]; end
%---
[ok,isrun] = soundmexpro('started');
if isrun, status(2) = 1; end
%---
if status(1) ~= status(2)
    soundmexpro('exit');
    status = [];
    error('status inconsistency - soundmexpro is stopped')
end
%---

ASIODriver = 'ASIO Fireface';%'ASIO Fireface USB';
SampFreq   = 44100;
LPC        = length(PlayChannelNums);
Tout       = 5;%sec 


if size(PlayMatrix,2) ~= LPC
    error('Invalid size of parameter PlayMatrix. PlayMatrix needs to have %d columns.',LPC);
end

switch Mode
    %----------------------------
    case {'init','Init'}
        if nargin < 2
            trigmode = 'Idx';
        end        
        if nargin > 2
            warning('MATLAB:EEGTrigPlay:TooManyArgInModeInit', 'arguments after ''init'' and ''trigmode'' are ignored' );
        end

        %---------------------------------------
        % initialize and start soundmexpro 
        %---------------------------------------
        ok = soundmexpro('exit');
        %----
        ok = soundmexpro('init', ...     % command name
            'force',      1, ...
            'driver',     ASIODriver, ...
            'output',     PlayChannelNums(:)');
        
        if ~ok, error('MATLAB:MRTPLAYREC:CannotInitializeSoundmex','error calling ''init'''); end
        %----
        ok = soundmexpro('cleardata');
        if ~ok, error('MATLAB:MRTPLAYREC:CannotCleardata','error calling ''cleardata'''); end
        %----
        status = 'initialized';
    case {'start','Start'}
        % starts playback without stopping
        ok = soundmexpro('start','length',0);
        if ~ok,  error(['error calling ''start''' ]); end
        %----
        % überprüfen ob/warten bis Wiedergabe+Aufnahme läuft
        isrun = 0; tout  = 0;
        tic
        while ~isrun
            if tout > Tout
                error('time out when starting soundmexpro for playback');
            end
            [ok,isrun] = soundmexpro('started');
            drawnow;
            tout = toc;
        end
        [ok, tload] = soundmexpro('trackload');
        if tload
            status = 'running';
        else
            status = 'waiting';
        end
    %----------------------------    
    case {'load','Load','data','Data'}
        if nargin < 2
           error('Parameter Playmatrix is needed in Mode %s',Mode); 
        end
        
        [ok, tload] = soundmexpro('trackload');
        if tload < 100
            status = 'loaded';
        else
            status = 'waiting';
            warning('wait till more blocks are already played');
            
        end
              
        if ~iscell(PlayMatrix)
            PlayMatrix = {PlayMatrix};
        end
        NumOfSignals = 0;
        for c = 1:length(PlayMatrix)
            for z = 1:size(PlayMatrix{c},3)
                NumOfSignals = NumOfSignals+1;
                
                if isempty(TrigID)
                    Trig.ID(NumOfSignals,1) = NumOfSignals;
                    Trig.IDMode             = 'Int';
                end
                [TrigInfo] = EEGTrigID2info(Trig.ID(NumOfSignals,:),Trig.IDMode,Trig.BitRef,Trig.NumOfBits,Trig.TwoBitScheme);
                
                if ischar(PlayMatrix{c})
                    PlayData = 0;
                else
                    PlayData{NumOfSignals}(:,1:2) = squeeze(PlayMatrix{c}(:,1:2,z));
                end
                TrigInfo.SigSamp = size(PlayData{NumOfSignals},1);
                TrigInfo.LatSamp = TrigLatSamp;
                TrigInfo.DurSamp = TrigDurSamp;

                PlayData{NumOfSignals}(:,3)   = GetTrigger(TrigInfo);
                                    
            end
        end
        
        %=======================
        for L = 1:BlockLoop
            %=======================
            if BlockOrder == -1
                NOS = randperm(NumOfSignals);
            elseif BlockOrder == 1
                NOS = 1:NumOfSignals;
            elseif length(BlockOrder) == NumOfSignals
                NOS = BlockOrder;
            else
                error('invalid value for BlockOrder');
            end
            TrigIDPlayed = Trig.ID(NOS);
            %=======================
            for nos = NOS
                ok = soundmexpro('loadmem', ...      % command name
                    'data', PlayData{nos}, ...          % data vector
                    'loopcount', SigLoop);
                if  ~ok, error('error calling ''loadmem'''); end
            
                JitterSig = zeros(randperm(JitterInBlockSamp,1),3);
                ok = soundmexpro('loadmem', ...      % command name
                    'data', JitterSig);
                if  ~ok, error('error calling ''loadmem'''); end                
            end
                        
            JitterSig = zeros(randperm(JitterOutBlockSamp,1),3);
            ok = soundmexpro('loadmem', ...      % command name
                'data', JitterSig);
            if  ~ok, error('error calling ''loadmem'''); end

            %=======================
        end %L
        
        for z = 1:size(PlayMatrix,3)
            ok = soundmexpro('loadmem', ...      % command name
                'data', PlayTrigMatrix, ...          % data vector
                'loopcount', 1);
            if  ~ok, error('error calling ''loadmem'''); end
        end
        
        % EEGPlayMatrix, TrigID, TrigT, NumOfLoops
        [ok,isrun] = soundmexpro('started');
        if isrun
            [ok,XRun] = soundmexpro('xrun');
            if XRun
                soundmexpro('exit');
                warning('MATLAB:PlayRecordLoopedASIO:DropOuts','Dropouts during playback/recording and irregular stop of soundmexpro')
            end
            
            %[ok,XRun] = soundmexpro('clipcount');
            
        else
            error('sounddevice is not running - please use Mode ''init'' before ''play''');
        end
    %----------------------------    
    case {'stop','Stop','exit'}
        if nargin >= 2
            warning('MATLAB:EEGTrigPlay:TooManyArgInModeStop', 'arguments after ''init'' are ignored' );
        end

        [ok] = soundmexpro('stop');
        pause(0.5);
        [ok] = soundmexpro('exit');
        
end

%========================================================================
function [Trigger] = GetTrigger(TI)

Trigger = zeros(TI.SigSamp,1);
Trigger(TI.LatSamp:TI.DurSamp) = TI.Trigword;



