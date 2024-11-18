function EEGTrigForwardUserResponse(answer)

%========================================================
%Forwarding response as Trigger to the EEG recording
%========================================================
% to be improved in future
% [Trig] =  EEGTrigID2info(0, [], eeg.System);
% eeg.NoB = Trig.NumOfBits;
global eeg

if isnumeric(answer)
    answer = num2str(answer);
end

%eeg.ForwResp2Trig
if  ~isnan(str2double(answer)) 
    RespTrigID = str2double(answer);
    
    %sending a 100 samples long Trigger to the EEG via RME->TriggerBox representing the given
    %user answer
    [Trig] =  EEGTrigID2info([0 RespTrigID], eeg.BitSplitScheme, eeg.System);
    disp(sprintf('UserResponse TriggerID: %d   <=>   [%d,%d]   <=>   %s',Trig.TrigInt, Trig.TrigIntMulti(1),Trig.TrigIntMulti(2),num2str(Trig.TrigBitmask)) );
    disp('-------------------------------------------------------------------------' );
    soundmexpro('loadmem', ...      % command name
        'data', [zeros(480,2),  Trig.TrigWord .*ones(480,1)] , ...  % data matrix
        'name', 'UserResp', ...            % name used in track view for this vector
        'loopcount', 1  );
end
