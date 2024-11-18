function [ trigstr, trigword, trigmask ] = TrigID2Trig( TrigID, TrigLED )
%
%
%
%
%  Copyright (C) 2017 Manfred Mauermann
%


if nargin < 2, TrigLED = []; end
if nargin < 1, help(mfilename); end

if ischar(TrigID)
    if TrigID(1)=='T'
        TrigID = str2double(TrigID(2:end));
        if TrigID == round(TrigID) & TrigID>=1 & TrigID<=12
            TrigID = 2^(TrigID-1);
        else
            error('TrigID needs to be an Integer between 1 and 4096 or a string ''T1'' ... ''T12'', or a cellarray with integers between 1 and 12');
        end
    end
end
if iscell(TrigID)   
        TrigID = unique([TrigID{:}]);
        if any(TrigID == round(TrigID)) & min(TrigID)>=1 & max(TrigID)<=12
            TrigID = sum(2.^([TrigID]-1));
        else
            error('TrigID needs to be an Integer between 1 and 4096 or a string ''T1'' ... ''T12'', or a cellarray with integers between 1 and 12');
        end
    
end

TrigID = round(TrigID);
if TrigID > 4095 || TrigID < 1
    error('TrigID needs to be an Integer between 1 and 4096 or a string ''T1'' ... ''T12''');
end


%Triggerwort fuer Triggerbox-NEW
% zeros according to the size of the needed bitmask
trigmask  = zeros(1,32);

tmp = de2bi(TrigID);
trigmask(1:length(tmp)) = trigmask(1:length(tmp))+tmp;

if ~isempty(TrigLED)
    trigmask(21:end) = trigmask(21:end) + TrigLED;
end

trigstr  = strrep(num2str(trigmask),' ','');
trigword = bin2dec(trigstr)/(1*2147483648);
%trigword = bin2dec(trigstr)/(1*16777215);
