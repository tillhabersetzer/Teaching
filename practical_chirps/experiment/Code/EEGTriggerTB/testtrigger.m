SF = 44100
soundmexpro('exit')
if 1 ~= soundmexpro('init', ...
        'driver', 1, ...    % driver index
        'output', [8], ...    % list of output channels to use
        'ramplen', 256, ...
        'track', 1 ...
        )
    error(['error calling ''init''' error_loc(dbstack)]);
end


if 1 ~= soundmexpro('trackmap', ... % command name 
        'track', [0] ... % track map [sig_left sig_right noise_left noise_right trigger]
        ) 
    error(['error calling ''trackmap''' error_loc(dbstack)]);
end


% general Trigger related properties
RMEDELAY       = 43;% in samples, 43 the for RME ADI-8 DS due to the reconstruction filter
USERDELAY      = 0; %in samples
amplifierdelay = 0.044;% in ms, delay of one-sample-click response relative to trigger,
acousticdelay  = 1; %in ms

TrigInfo.Delay        = acousticdelay + amplifierdelay + (RMEDELAY + USERDELAY)*1000/SF;    % in ms
TrigInfo.SF           = SF;
TrigInfo.TrigDur      = 100;%ms





%[trigdata] = get_trigger(TrigID, TrigInfo)

if 1 ~= soundmexpro('start','length', 0), error(['error calling ''start''' error_loc(dbstack)]); end

TrigID = [1 12];% {[250,8]};%[1 3 5 7 ];%[1 3 5 7 9] ;% {[111,0]};% [1 3 5 7 9] 

[TrigInfo] = EEGTrigID2info(TrigID)
% TrigInfo.Trigword = 1
%TrigWord = sum((2.^(32-[12]))/(2^32)); sum((2.^(32-[12]))/(2^32));
%TrigWord = 1;% 
%TrigWord= sum((2.^(-[12 8])))
data = ones(2048,1).*TrigInfo.Trigword;
% 
amplitude = bin2dec(    '1000 0000    0001  0000     0000 0000');
amplitude = ((amplitude+1)/(16777216))
data = ones(2048,1).*amplitude


[success] =  soundmexpro('loadmem', 'data', data, ...
    'track',0, ...
    'loopcount', 80); % play single loop
if success ~= 1, error(['error calling ''loadmem''' error_loc(dbstack)]);   end
 