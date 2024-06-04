function [TrigInfo] = EEGTrigID2info_(TrigID,BitSplitScheme,System,TrigDur,ASIO,BitRef,NumOfBits)
%
%[TrigInfo] =  EEGTrigID2info(TrigID,BitSplitScheme,System[,PlayFlag,BitRef,NumOfBits])
%
% EEGTrigID2info(TrigID) translates an EEG Trigger ID for the
% actual triggerbox to other representations like the triggerword bitmask
% etc.
%
% Input:
% TrigID         - Trigger identifier (see below)
% BitSplitScheme - vector with x numbers of bits to split the complete bitmask
%                  with 'NumOfBits' bits into x sub-bitmasks. This might
%                  be used for example to split the 12 bits of the new
%                  trigger box into 8 bits which are only used for trigger
%                  labeling and 4 bits which are used as well to mark the
%                  LEDs on the response box - details see below. Default is
%                  [NumOfBits 0]
% System         - 'BioSemi' or 'SynAmps' is needed to consider differences in the
%                  Trigger value recognized by the respective recording software
%                  according to the given Triggerword - there appears to be a bitflip in the
%                  current SynAmps setup (low Triggerwords => high Trigger
%                  integer values).
%                  (default: 'Biosemi')
%
%-- optional parameters for testing purposes
% TrigDur        - if > 0 a trigger test mode will be started and the triggers 
%                  can be played with a duration of 'TrigDur' ms 
%                  (default: 0)
% ASIO           - structure with the fields: Driver, SampFreg, TrigChan
%                  needed  in the test mode to access (via soundmexpro) the used 
%                  sound device 'Driver' with
%                  samplingfrequency 'SampFreq' via the 
%                  ASIO (SPDIF) channel 'TrigChan'    
%                  (default: ASIO.Driver   = 'ASIO Fireface USB';
%                            ASIO.SampFreq = 48000;
%                            ASIO.TrigChan = 8;) 
% -- optional profi parameters - use them only when you know what you are doing ---
% BitRef         - maximum number of bits used by the specific triggerbox
% NumOfBits      - number of bits used for the specific EEG-System trigger
%                  (default: Biosemi - 12 bit, Synamps - 8 bit)
%
% There are three modes to describe the Trigger IDs 'TrigID': 'Bitmask',
% 'LED-numbers' and 'Integer' depending on the data format used for
% TrigID:
% 1) LED-numbers - to generate a triggerword explicitly lighting specific LEDs
%    on the new trigger box. Use a cell array which contains a vector with
%    the LEDs that should light up for the trigger word to be generated.
%    (e.g. {[1 3 5 12]} will generate a trigger word for which the LEDs with the numbers
%     1, 3, 5, and 12 will light up). If you want to generate multiple
%     triggerwords at once in this mode use a column cell array of row vectors
%     (e.g. {[1 12];[1 2 3]; ...} ). This mode might be used if an online
%     control of the measurement via the flashing LEDs should be easily traceable.
% 2) Bitmask - to describe the wanted trigger via an bitmask give the TrigID
%    as respective binary number(s) noted as string or for multiple triggers
%    as character array. The number of digits is given by the parameter NumOfBits.
%    Examples: TrigID = '100000000010' will generate a 12 bit
%    Triggerword that would light up the LEDs 1 and 11, and is related to a
%    integer value of 2050 (2048 + 2). Multiple trigger ID can be given as
%    character array where different triggers are given in different rows
%    like TrigID = ['100000000010';'000000000011';'110000000000'] will
%    generate three triggers light up the LEDs 1 and 11, 11 and 12 or, 1 and
%    2 respectivley and are related to the integer values of 2050, 3, and
%    3072.
% 3) Integer - by giving a single integer number between 1 and
%    2^NumOfBits-1, or for multiple triggers as column vector of integers,
%    will set the rspective trigger bits according to the integer value.
%    Example: TrigID = 2050 would light up the LEDs 1 and 11 and is related
%    to a bitmask '100000000010'
%
% The 'BitSplitScheme' allows to divide the 'NumOfBits' in subsets of
% bitranges. This might be used to access via the triggers directly
% different groups of stimuli with a several stimuli in it- Or to separate
% between bits accessing the response box LEDs and bits used for triggering
% alone. The default is [NumOfBits 0]; other values might be [8 4] (giving
% the number of bits of the low "byte" to the high "byte" from left to
% right) to separate the 12 bits of the new trigger box in ranges covering trigger
% alone (8bit) and response box LED control (4 bit). The vectors elements
% need to sum up to NumOfBits. The TrigIDs are separated in respective
% subsets. In case of NumOfBits = 12 and BitSplitScheme = [8 4] for example
% TrigID (as LED-numbers) = {[1 8],[1 4];[2 7],[2 3]} would generate the two
% triggers with the respective splitted TrigID integer scheme [129 9;66 6] given
% by TrigInfo.TrigIntMulti and vice versa (see TrigInfo.TrigLEDMulti).
% The integer scheme for the whole bitrange (in ths case 12 bit)
% is given in TrigInfo.TrigInt as [2073;1062]. Bitmask and Trigword are
% identical for splitted and not splitted bit scheme.
%
%
%
% (c) Manfred Mauermann 2017
%



% 20.7.2018 include Triggercheck via Playback over ASIO device
%
%

% 09.02.2018 included MEG and make the BitRef default system dependent
%The BitRef value might be used to create small offsets in the TriggerWord
%changed order of parameters
%



% check arguments
if nargin < 2, BitSplitScheme = []; end
if nargin < 3, System         = []; end
if nargin < 4, TrigDur        = []; end
if nargin < 5, ASIO           = ''; end
if nargin < 6, BitRef         = []; end
if nargin < 7, NumOfBits      = []; end

%if nargin < 2, TrigIDType    = []; end
if nargin < 1, error('Parameter TrigID is needed') ; end

if ischar(TrigID)
    TrigIDType = 'bitmask';
elseif iscell(TrigID)
    TrigIDType = 'led';
elseif isnumeric(TrigID)
    TrigIDType = 'int';
else
    error('Invalid TrigIDType');
end

% set default values if needed
if isempty(System),            System     = 'BioSemi';  end %SynAmps, MEG, BioSemi
if isempty(TrigDur),           TrigDur    = 0; end
if isempty(ASIO),           ASIO.Driver   = 'ASIO Fireface'; 
                            ASIO.SampFreq = 48000;
                            ASIO.TrigChan = 12; %ASIO.TrigChan = 8; 
end


switch System
    case {'MEG','meg','triux','standard','default'}
        if isempty(BitRef),         BitRef         = 24;     end
        if isempty(NumOfBits),      NumOfBits      = 12;     end
        if isempty(BitSplitScheme), BitSplitScheme = [NumOfBits]; end
        SetLowBits = false;
    case {'BioSemi','biosemi'}
        if isempty(BitRef),         BitRef         = 24;     end
        if isempty(NumOfBits),      NumOfBits      = 12;     end
        if isempty(BitSplitScheme), BitSplitScheme = [NumOfBits]; end
        %BitSplitScheme = fliplr(BitSplitScheme);
        SetLowBits = false;
    case {'BioSemi_OldTB16','biosemi_OldTB16'}
        if isempty(BitRef),         BitRef         = 24;     end
        if isempty(NumOfBits),      NumOfBits      = 16;     end
        if isempty(BitSplitScheme), BitSplitScheme = [NumOfBits]; end
        %BitSplitScheme = fliplr(BitSplitScheme);
        SetLowBits = false;
       
    case {'Synamps'}
        if isempty(BitRef),         BitRef         = 16;     end
        if isempty(NumOfBits),      NumOfBits      = 8;      end
        if isempty(BitSplitScheme), BitSplitScheme = [NumOfBits]; end
        SetLowBits = false;
    otherwise
        error('selected ''System'' is not available for EEGTrigID2info')
end
LowBitMask = logical(ones(1,BitRef-NumOfBits));
if sum(BitSplitScheme) ~= NumOfBits | size(BitSplitScheme,2) > NumOfBits | size(BitSplitScheme,2) < 1
    error('BitSplitScheme has to be a 1 up to ''NumOfBits'' element row vector summing up to ''NumOfBits'' (here %d)',NumOfBits)
end

switch TrigIDType
    case {'LED','led','idx','Idx','IDX','Index'}
        if size(TrigID,2) ~= size(BitSplitScheme,2)
            error('size of 2nd dimension of the TrigID cell array (%d) needs to match with the length of BitSplitScheme (%d)',size(TrigID,2),size(BitSplitScheme,2));
        end
        
        for t = 1:size(TrigID,1)
            BM_  = logical([]);
            for g = 1:size(TrigID,2)
                bm = false( 1,BitSplitScheme(g) );
                trigid = unique(TrigID{t,g});
                if any(trigid > BitSplitScheme(g)) | any(trigid < 0)
                    error('Invalid LED indices in TrigID(%d,%d) values in index vector ([%s]) has to be in the range 1-%d',t,g,sprintf('%d ',TrigID{t,g}),BitSplitScheme(g));
                end
                bm(trigid) = true;
                BM_        = [BM_ bm];
            end
            TrigID_(t,:) = strrep(num2str(logical(BM_)),' ','');
        end
        
    case {'bit','Bit','bitmask','Bitmask','BitMask'}
        if  size(TrigID,2) ~= NumOfBits
            error('Length of bitmask (%d) differ from nominal NumOfBits (%d) for System ''%s'' ',size(TrigID,2),NumOfBits,System);
        end
        
        TrigID_  = TrigID;
        
    case {'int','Int','integer','Integer'}
        
        if size(TrigID,2) ~= size(BitSplitScheme,2)
            error('size of 2nd dimension of the TrigID integer array (%d) needs to match with the length of BitSplitScheme (%d)',size(TrigID,2),size(BitSplitScheme,2));
        end
        for t = 1:size(TrigID,1)
            BM_  = '';
            for g = 1:size(BitSplitScheme,2)
                if any(TrigID(:,g) > 2.^(BitSplitScheme(g))-1) | any(TrigID(:,g) < 0)
                    error('Values for TrigID in mode %s(%d) has to be between 0 and %d',TrigIDType,BitSplitScheme(g), 2^(BitSplitScheme(g))-1);
                end
                BM_    = [dec2bin(TrigID(t,g),BitSplitScheme(g)) BM_];
                
            end
            TrigID_(t,:) = BM_;
        end
end

%=======================================================================
switch System
    %---------------------------------------------
    case {'MEG','meg','triux','standard','default'}
        Bitmask = logical(TrigID_-'0');
        TrigInfo.TrigBitmask  = Bitmask;
        [x,y]                 = find(Bitmask);
        for k = 1:max(x)
            TrigInfo.TrigLED{k,1} = uint8(y(x==k))';
        end
        TrigInfo.TrigInt      = uint16(bin2dec(num2str(Bitmask)));
        TrigInfo.System       = System;
        
        if SetLowBits
            TrigInfo.TrigWord  = sum( bitset(0, repmat(BitRef:-1:1,size(Bitmask,1),1), [Bitmask LowBitMask])./(2.^BitRef-1),2);
        else
            TrigInfo.TrigWord  = sum( bitset(0, repmat(BitRef:-1:BitRef-NumOfBits+1,size(Bitmask,1),1), Bitmask )./(2.^BitRef-1),2);
        end
        for t = 1:size(Bitmask,1)
            for g = 1:size(BitSplitScheme,2)
                TrigInfo.TrigLEDMulti{t,g}  = find( Bitmask(t, sum(BitSplitScheme(1:g-1))+1:sum(BitSplitScheme(1:g)) ) );
                TrigInfo.TrigIntMulti(t,g)  = bin2dec( strrep( num2str(  Bitmask(t, sum(BitSplitScheme(1:g-1))+1:sum(BitSplitScheme(1:g)) ) ),' ','') );
            end
        end
        %---------------------------------------------
    case {'BioSemi','biosemi','Biosemi'}
        Bitmask = logical(TrigID_-'0');
        switch TrigIDType
            case {'LED','led','idx','Idx','IDX','Index'}
                Bitmask = fliplr(Bitmask);
        end
        TrigInfo.TrigBitmask   = Bitmask;%as they appear in the BioSemi system
        TrigInfo.TrigInt       = uint16(bin2dec(num2str(Bitmask)));%int values as they appear in the BioSemi system
        TrigInfo.System        = System;
        TrigInfo.DeviceBitmask = logical(dec2bin(TrigInfo.TrigInt,16)-'0');
        
        %Bitmask as it appears in the Biosemi system
        for t = 1:size(TrigInfo.DeviceBitmask,1)
            TrigInfo.DeviceInt(t,:)  = [ uint8(bin2dec(num2str(TrigInfo.DeviceBitmask(t,9:16)))), uint8(bin2dec(num2str(TrigInfo.DeviceBitmask(t,1:8))))];
        end
        
        [x,y]                  = find(fliplr(Bitmask));
        for k = 1:max(x)
            TrigInfo.TrigLED{k,1} = uint8(y(x==k))';
        end
        
        LEDBitmask = circshift(TrigInfo.DeviceBitmask,-4,2);
        LEDBitmask = LEDBitmask(:,1:12);
        for t = 1:size(Bitmask,1)
            for g = 1:size(BitSplitScheme,2)
                TrigInfo.TrigLEDMulti{t,g}  = find( LEDBitmask (t, 13-(sum(BitSplitScheme(1:g-1))+1:sum(BitSplitScheme(1:g))) ) );
                TrigInfo.TrigIntMulti(t,g)  = bin2dec( num2str( fliplr( LEDBitmask (t, 13-(sum(BitSplitScheme(1:g-1))+1:sum(BitSplitScheme(1:g))) ) )  ) ) ;
            end
        end
        
        LowBitMaskWord     = bin2dec('00000000 0000 0000 0000 1111')/(2.^BitRef-1);%bin2dec('00000000 0000 0000 0000 1111')/(1*16777215)
        TrigInfo.TrigWord  = ( sum( bitset(0, repmat(BitRef-NumOfBits+1:1:BitRef,size(Bitmask,1),1), Bitmask )./(2.^BitRef-1),2) ...
            + any(Bitmask(:,1:4),2).* LowBitMaskWord);
        %---------------------------------------------
    case {'BioSemi_OldTB16','biosemi_OldTB16'}
        Bitmask = logical(TrigID_-'0');
        switch TrigIDType
            case {'LED','led','idx','Idx','IDX','Index'}
                Bitmask = fliplr(Bitmask);
        end
        TrigInfo.TrigBitmask   = Bitmask;%as they appear in the BioSemi system
        TrigInfo.TrigInt       = uint16(bin2dec(num2str(Bitmask)));%int values as they appear in the BioSemi system
        TrigInfo.System        = System;
        TrigInfo.DeviceBitmask = logical(dec2bin(TrigInfo.TrigInt,16)-'0');
        
        %Bitmask as it appears in the Biosemi system
        for t = 1:size(TrigInfo.DeviceBitmask,1)
            TrigInfo.DeviceInt(t,:)  = [ uint8(bin2dec(num2str(TrigInfo.DeviceBitmask(t,9:16)))), uint8(bin2dec(num2str(TrigInfo.DeviceBitmask(t,1:8))))];
        end
        
        [x,y]                  = find(fliplr(Bitmask));
        for k = 1:max(x)
            TrigInfo.TrigLED{k,1} = uint8(y(x==k))';
        end
        
        for t = 1:size(Bitmask,1)
            for g = 1:size(BitSplitScheme,2)
               TrigInfo.TrigIntMulti(t,g)  = bin2dec( strrep( num2str(  Bitmask(t, sum(BitSplitScheme(1:g-1))+1:sum(BitSplitScheme(1:g)) ) ),' ','') );
            end
        end
        
        TrigInfo.TrigWord  =  sum( bitset(0, repmat(NumOfBits:-1:1,size(Bitmask,1),1), Bitmask )./(2.^(BitRef-1)),2) ;
        
        
        %---------------------------------------------
    case {'Synamps','synamps','Syn','syn','SynAmps','synAmps'} % for any reason (most probably cable) the bitconfiguration for the synamps ist flipped. To keep it consistent only the Triggerword is adapted respectively
        TrigInfo.TrigBitmask  = Bitmask;
        [x,y]                 = find(Bitmask);
        for k = 1:max(x)
            TrigInfo.TrigLED{k,1} = uint8(y(x==k))';
        end
        TrigInfo.TrigInt      = uint16(bin2dec(num2str(Bitmask)));
        TrigInfo.System       = System;
        
        TrigInfo.TrigWord = sum( bitset(0, repmat(BitRef-NumOfBits+1:BitRef,size(Bitmask,1),1), Bitmask )./(2.^BitRef-1),2);
        for t = 1:size(Bitmask,1)
            for g = 1:size(BitSplitScheme,2)
                TrigInfo.TrigLEDMulti{t,g}  = find( Bitmask(t, sum(BitSplitScheme(1:g-1))+1:sum(BitSplitScheme(1:g)) ) );
                TrigInfo.TrigIntMulti(t,g)  = bin2dec( strrep( num2str(  Bitmask(t, sum(BitSplitScheme(1:g-1))+1:sum(BitSplitScheme(1:g)) ) ),' ','') );
            end
        end
        
    otherwise
        error('selected ''System'' is not available for EEGTrigID2info')
end


try
    
    NumOfTW = size(TrigInfo.TrigWord,1);
    
    
    if TrigDur > 0
        SampFreq       = 48000;
        TrigDuration   = 100;%ms
        TriggerChannel = 8;
        
        
        TrigSig                    = ones(TrigDur/1000*SampFreq,1) * TrigInfo.TrigWord';
        
        ok = soundmexpro('init', ...     % command name
            'force',      1, ...
            'driver',     ASIO.Driver, ...
            'samplerate', ASIO.SampFreq, ...
            'output',     ASIO.TrigChan, ...
            'track',      1, ...
            'input',      [], ...
            'numbufs',    1, ...
            'autocleardata', 1, ...
            'ramplen', 0);
        
        ok = soundmexpro('loadmem', ...      % command name
            'data', zeros(48,1), ...          % data vector
            'loopcount', 1 ...
            );
        
        ok = soundmexpro('start','length',0);
        if ~ok
            error(['error calling ''start''' ]);
        end
        
        TrigNum  = 0;
        
        tn_str =  sprintf('%d, ',1:NumOfTW);
        pr_str = sprintf('\ntype a number an then press RETURN:\n %s to play the respective Trigger from the list, \n 0 to toggle or \n-1 to exit\nthe trigger test playback \n', tn_str);
        while TrigDur > 0
            k = input(pr_str);
            switch k
                case -1
                    TrigDur  = 0;
                    TrigNum  = 0;
                case 0
                    TrigNum = TrigNum + 1;
                    TrigNum = mod(TrigNum-1, NumOfTW) + 1;
                    
                case num2cell(1:NumOfTW)
                    TrigNum = k;
            end

            if TrigNum
                disp(sprintf('playing trigger number %d: %s',TrigNum, num2str(TrigInfo.TrigBitmask(TrigNum,:))));
                ok = soundmexpro('loadmem', ...      % command name
                    'data', TrigSig(:,TrigNum), ...          % data vector
                    'loopcount', 1 ...
                    );
                if ~ok
                    error(['error calling ''loadmem''' ]);
                end
                pause(0.01);
            end            
        end %while
        
    end %if TrigDur > 0
catch
    soundmexpro('exit');
end
soundmexpro('exit');

% switch TRIG
%     case 1
%         whichtrig = bin2dec('10000000 0000 0000 0000 0000')/(1*16777215);   %TRIG1
%     case 2
%         whichtrig = bin2dec('01000000 0000 0000 0000 0000')/(1*16777215);   %TRIG2
%     case 3
%         whichtrig = bin2dec('00100000 0000 0000 0000 0000')/(1*16777215);   %TRIG3
%     case 4
%         whichtrig = bin2dec('00010000 0000 0000 0000 0000')/(1*16777215);   %TRIG4
%     case 5
%         whichtrig = bin2dec('00001000 0000 0000 0000 0000')/(1*16777215);   %TRIG5
%     case 6
%         whichtrig = bin2dec('00000100 0000 0000 0000 0000')/(1*16777215);   %TRIG6
%     case 7
%         whichtrig = bin2dec('00000010 0000 0000 0000 0000')/(1*16777215);   %TRIG7
%     case 8
%         whichtrig = bin2dec('00000001 0000 0000 0000 0000')/(1*16777215);   %TRIG8
%     case 9
%         whichtrig = bin2dec('00000000 1000 0000 0000 1111')/(1*16777215);   %LED1 + TRIG9
%     case 10
%         whichtrig = bin2dec('00000000 0100 0000 0000 1111')/(1*16777215);   %LED2 + TRIG10
%     case 11
%         whichtrig = bin2dec('00000000 0010 0000 0000 1111')/(1*16777215);   %LED3 + TRIG11
%     case 12
%         whichtrig = bin2dec('00000000 0001 0000 0000 1111')/(1*16777215);   %LED4 + TRIG12
% end