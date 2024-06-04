classdef EEGTrigger
    
    properties (SetAccess = protected)
        BitSplitScheme;
        NumOfBits;
        BitRef;
        System;
        ID_LED;
        ID_Int;
        ID_BitMask;
        ID_BitStr;
        ID_LEDMulti;
        ID_IntMulti;
        Word
    end
    properties (Access = public, Hidden = true)
        TrigID
    end
    
    methods
        %% ========================================================================
        function trig = EEGTrigger(TrigID,varargin)
            %parse input of constructor
            p = inputParser;
            default.System         = 'Biosemi';
            default.BitSplitScheme = 12;
            default.BitRef         = 24;
            default.NumOfBits      = 12;
            
            expected.Systems       = {'BioSemi','biosemi','Biosemi','standard','default', ...
                'Synamps','synamps','Syn','syn','SynAmps','synAmps'};
            
            addParameter(p,'System',        default.System,        @(x) any(validatestring(x,expected.Systems)));
            addParameter(p,'BitSplitScheme',default.BitSplitScheme,@isnumeric);
            addParameter(p,'BitRef',        default.BitRef,        @isnumeric);
            addParameter(p,'NumOfBits',     default.NumOfBits,     @isnumeric);
            addRequired( p,'SetTrigID');
            parse(p,TrigID,varargin{:});
            
            fn = fieldnames(p.Results);
            fn(strmatch('SetTrigID',fn)) = [];
            for k =  1:length(fn), trig.(fn{k}) = p.Results.(fn{k}); end
            
            % additional checks
            if sum(trig.BitSplitScheme) ~= trig.NumOfBits | size(trig.BitSplitScheme,2) > trig.NumOfBits | size(trig.BitSplitScheme,2) < 1
                error('BitSplitScheme has to be a 1 up to ''NumOfBits'' element row vector summing up to ''NumOfBits'' (here %d)',trig.NumOfBits)
            end
            %trig = SetTrigID(trig,TrigID);
            trig.TrigID = TrigID;
        end
        
        function trig = set.TrigID(trig,TrigID) %trig = SetTrigID(trig,TrigID)
            if  iscell(TrigID) %LED-Index
                if size(TrigID,2) ~= size(trig.BitSplitScheme,2)
                    error('size of 2nd dimension of the TrigID cell array (%d) needs to match with the length of BitSplitScheme (%d)',size(TrigID,2),size(trig.BitSplitScheme,2));
                end
                trig.TrigID = TrigID;
                Bitmask = LED2bitmask(trig,TrigID);
                trig    = bitmask2all(trig,Bitmask);
            elseif ischar(TrigID) %BitStr
                if  size(TrigID,2) ~= trig.NumOfBits
                    error('Length of bitmask (%d) differ from nominal trig.NumOfBits (%d) for System ''%s'' ',size(TrigID,2),trig.NumOfBits,trig.System);
                end
                trig.TrigID = TrigID;
                Bitmask = logical(TrigID-'0');
                trig    = bitmask2all(trig,Bitmask);
            elseif islogical(TrigID)
                if size(TrigID) ~= size(trig.BitSplitScheme(trig.BitSplitScheme~=0),2)
                    error('size of 2nd dimension of the TrigID integer array (%d) needs to match with the length of BitSplitScheme (%d)',size(TrigID,2),size(trig.BitSplitScheme,2));
                end
                trig.TrigID = TrigID;
                Bitmask = TrigID;
                trig    = bitmask2all(trig,Bitmask);
            elseif isnumeric(TrigID)
                if size(TrigID) ~= size(trig.BitSplitScheme(trig.BitSplitScheme~=0),2)
                    error('size of 2nd dimension of the TrigID integer array (%d) needs to match with the length of BitSplitScheme (%d)',size(TrigID,2),size(trig.BitSplitScheme,2));
                end
                trig.TrigID = TrigID;
                Bitmask = int2bitmask(trig,TrigID);
                trig    = bitmask2all(trig,Bitmask);
            end
            
            
        end
        
    end
    %%
    methods (Access = protected, Hidden = true)
        function trig = bitmask2all(trig,Bitmask)
            [x,y]           = find(Bitmask);
            for k = 1:max(x)
                ID_LED{k,1} = uint8(y(x==k))';
            end
            trig.ID_LED       = ID_LED;
            trig.ID_Int       = uint16(bin2dec(num2str(Bitmask)));
            trig.ID_BitMask   = Bitmask;
            trig.ID_BitStr    = '';
            for t = 1:size(Bitmask,1)
                trig.ID_BitStr(t,:)    = strrep( num2str(Bitmask(t,:)),' ','');
                for g = 1:size(trig.BitSplitScheme,2)
                    trig.ID_LEDMulti{t,g}  = find( Bitmask(t, sum(trig.BitSplitScheme(1:g-1))+1:sum(trig.BitSplitScheme(1:g)) ) );
                    trig.ID_IntMulti(t,g)  = bin2dec( strrep( num2str(  Bitmask(t, sum(trig.BitSplitScheme(1:g-1))+1:sum(trig.BitSplitScheme(1:g)) ) ),' ','') );
                end
            end
            
            switch trig.System
                case {'BioSemi','biosemi','Biosemi','standard','default'}
                    trig.Word = sum( bitset(0, repmat(trig.BitRef:-1:trig.BitRef-trig.NumOfBits+1,size(Bitmask,1),1), Bitmask )./(2.^trig.BitRef-1),2);
                case {'Synamps','synamps','Syn','syn','SynAmps','synAmps'} % for any reason (most probably cable) the bitconfiguration for the synamps ist flipped. To keep it consistent only the Triggerword is adapted respectively
                    trig.Word = sum( bitset(0, repmat(trig.BitRef-trig.NumOfBits+1:trig.BitRef,size(Bitmask,1),1), Bitmask )./(2.^trig.BitRef-1),2);
                otherwise
                    error('selected ''System'' is not available for EEGTrigID2info')
            end
            
        end
        %------------------------------------------
        function Bitmask = int2bitmask(trig,ID)
            if size(ID,2) ~= size(trig.BitSplitScheme(trig.BitSplitScheme~=0),2)
                error('size of 2nd dimension of the TrigID integer array (%d) needs to match with the length of BitSplitScheme (%d)',size(TrigID,2),size(BitSplitScheme,2));
            end
            for t = 1:size(ID,1)
                BM_  = '';
                for g = 1:size(trig.BitSplitScheme,2)
                    if any(ID(:,g) > 2.^(trig.BitSplitScheme(g))-1) | any(ID(:,g) < 0)
                        error('Values for TrigID in mode %s(%d) has to be between 0 and %d',trig.IDType,trig.BitSplitScheme(g), 2^(trig.BitSplitScheme(g)-1));
                    end
                    BM_    = [BM_ dec2bin(ID(t,g),trig.BitSplitScheme(g))];
                    
                end
                TrigID_(t,:) = BM_;
            end
            Bitmask = logical(TrigID_-'0');
        end
        %------------------------------------------
        function Bitmask = LED2bitmask(trig,ID)
            if size(ID,2) ~= size(trig.BitSplitScheme,2)
                error('size of 2nd dimension of the TrigID cell array (%d) needs to match with the length of BitSplitScheme (%d)',size(trig.TrigID,2),size(trig.BitSplitScheme,2));
            end
            
            for t = 1:size(ID,1)
                BM_  = logical([]);
                for g = 1:size(ID,2)
                    bm = false( 1,trig.BitSplitScheme(g) );
                    trigid = unique(ID{t,g});
                    if any(trigid > trig.BitSplitScheme(g)) | any(trigid < 0)
                        error('Invalid LED indices in TrigID(%d,%d) values in index vector ([%s]) has to be in the range 1-%d',t,g,sprintf('%d ',TrigID{t,g}),trig.BitSplitScheme(g));
                    end
                    bm(trigid) = true;
                    BM_        = [BM_ bm];
                end
                TrigID_(t,:) = strrep(num2str(logical(BM_)),' ','');
            end
            Bitmask = logical(TrigID_-'0');
        end
        
        %------------------------------------------
        
    end
    
end

