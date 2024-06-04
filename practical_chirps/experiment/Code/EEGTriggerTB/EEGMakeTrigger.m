function [Trigger] = EEGMakeTrigger(TrigID, TrigDuration, TrigDelay, SampFreq, Bitsplit)
%
%[Trigger] = EEGMakeTrigger(TrigID, TrigDuration, TrigDelay, SampFreq)
%
%
% (c) Manfred Mauermann 2017

[Trigger] = EEGTrigID2info(TrigID,Bitsplit,'biosemi')
Trigger.TrigLED{1}

Trigger.DelaySamp = round(SampFreq/1000.* TrigDelay );
Trigger.LenSamp   = round(SampFreq/1000.* TrigDuration ); % in samples
Trigger.Duration  = Trigger.LenSamp./SampFreq;

Trigger.Sig        = zeros(Trigger.LenSamp+Trigger.DelaySamp,1);
Trigger.Sig(Trigger.DelaySamp+1:Trigger.DelaySamp+Trigger.LenSamp,1) ... 
                   = Trigger.TrigWord;



        