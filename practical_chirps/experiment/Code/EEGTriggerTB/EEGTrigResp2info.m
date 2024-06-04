function [RespInfo] = EEGTrigResp2info(data,System,BitSplitScheme,NumOfBits,BitRef)
%
% [RespInfo] = EEGTrigResp2info(input,System,BitSplitScheme
%
% Same as EEGTrigVal2info except that it uses the triggerword 
% of the Bio-Semi as input(as maybe obtained from the recording data).
% which is the integer trigger value devided by 2^32;
%
% see also: EEGTrigVal2info, EEGTrigID2info, int2bitmask
%
% (c) Manfred Mauermann 2018
%

 
if nargin < 5, BitRef         = [];    end
if nargin < 4, NumOfBits      = []; end
if nargin < 3, BitSplitScheme = []; end
 
switch System
%     case {'MEG','meg','triux','standard','default'}
%         if isempty(BitRef),         BitRef         = 24;     end
%         if isempty(NumOfBits),      NumOfBits      = 12;     end
%         if isempty(BitSplitScheme), BitSplitScheme = [NumOfBits]; end
%         SetLowBits = false;
    case {'BioSemi','biosemi'}
        if isempty(BitRef),         BitRef         = 24;     end
        if isempty(NumOfBits),      NumOfBits      = 12;     end
        if isempty(BitSplitScheme), BitSplitScheme = [NumOfBits]; end
        Bitmask                                    = dec2bin(data.*((2^BitRef-1)));
        TrigBitmask                                = Bitmask(1:NumOfBits);
        RespBitmask                                = Bitmask(NumOfBits+1:NumOfBits+4);
        [TrigInfo]                                 = EEGTrigID2info(TrigBitmask,[],'BioSemi');
        [RespInfo]                                 = EEGTrigID2info(strrep(num2str(TrigInfo.TrigBitmask),' ',''),BitSplitScheme,'BioSemi');
         RespInfo.RespBitmask                      = logical(RespBitmask-'0');   
         RespInfo.RespID                           = find(RespInfo.RespBitmask,4);  
%     case {'BioSemi_OldTB16','biosemi_OldTB16'}
%         if isempty(BitRef),         BitRef         = 24;     end
%         if isempty(NumOfBits),      NumOfBits      = 16;     end
%         if isempty(BitSplitScheme), BitSplitScheme = [NumOfBits]; end
%         TrigID                                     = round(val.*(2.^(NumOfBits-1)));
%         
        
    case {'Synamps'}
        if isempty(BitRef),         BitRef         = 24;     end
        if isempty(NumOfBits),      NumOfBits      = 16;      end
        if isempty(BitSplitScheme), BitSplitScheme = [NumOfBits]; end
        Bitmask                                    = dec2bin(data.*((2^BitRef)),16);
        TrigBitmask                                = fliplr(Bitmask(1:NumOfBits));
        [TrigInfo]                                 = EEGTrigID2info(TrigBitmask,[],'Synamps');
        [RespInfo]                                 = EEGTrigID2info(strrep(num2str(TrigInfo.TrigBitmask),' ',''),BitSplitScheme,'Synamps');
        RespBitmask                                = fliplr(Bitmask(NumOfBits+1+8:NumOfBits+4+8));
        RespInfo.RespBitmask                       = logical(RespBitmask-'0');   
        RespInfo.RespID                            = find(RespInfo.RespBitmask,4);  
    otherwise
        error('selected ''System'' is not available for EEGTrigID2info')
end
 
