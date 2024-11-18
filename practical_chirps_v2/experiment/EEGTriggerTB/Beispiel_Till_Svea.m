Liste = 3
Satz  = 5


ListeBits = 8;
SatzBits  = 8;

% if Liste > 2^ListeBits-1
%     error()
% end
% if Satz > 2^SatzBits;
%     error()
% end

[TrigInfo] = EEGTrigID2info([Liste Satz],[ListeBits SatzBits],'meg');

% bitmask_str = strrep(num2str(TrigInfo.TrigBitmask),' ','')
% [TrigInfo]  = EEGTrigID2info(bitmask_str(1:ListeBits+SatzBits),[ListeBits+SatzBits ],'meg')
% 
% 
% [TrigInfo] = EEGTrigID2info([TrigInfo.TrigIntNom],[ListeBits+SatzBits ],'meg')
% 
% 
% 
% 
%  strrep(num2str(TrigInfo.TrigBitmask),' ','')
