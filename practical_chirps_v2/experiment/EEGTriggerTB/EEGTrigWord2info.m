function [TrigInfo] = EEGTrigResp2info(input)
%
% [TrigInfo] = EEGTrigWord2info(trigword)
%
% Same as EEGTrigVal2info except that it uses the triggerword 
% of the Bio-Semi as input(as maybe obtained from the recording data).
% which is the integer trigger value devided by 2^32;
%
% see also: EEGTrigVal2info, EEGTrigID2info, int2bitmask
%
% (c) Manfred Mauermann 2017
%

% corrval = 1.9073486328125e-06;

corrval = 0;
Bitmask    = int2bitmask(round( (input-corrval).*2^32), 32);
TrigID     = find(Bitmask(1:12));
[TrigInfo] = EEGTrigID2info(TrigID);
TrigInfo.RespListID = find(Bitmask(16:-1:13));%Keys numbered from left to right as placed on the trigger box
