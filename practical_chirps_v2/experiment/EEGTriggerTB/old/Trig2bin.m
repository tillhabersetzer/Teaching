function [ trigstr, trigword, trigV ] = Trig2bin( Trig1, TrigLED )
%
%
%
%

if nargin < 2, TrigLED = []; end
if nargin < 1, help(mfilename); end

Trig1 = round(Trig1);
if Trig1 > 4096 || Trig1 < 1
    error('TriggerID needs to be an Integer between 1 and 4096');
end


%Triggerwort fuer Triggerbox-NEW
trigV  = zeros(1,24);

tmp = de2bi(Trig1);
trigV(1:length(tmp)) = trigV(1:length(tmp))+tmp;

if ~isempty(TrigLED)
    trigV(21:end) = trigV(21:end) + TrigLED;
end

trigstr  = strrep(num2str(trigV),' ','');
trigword = bin2dec(trigstr)/(1*16777215);