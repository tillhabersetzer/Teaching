function [answer] = EEGResponseBoxButton(TimeOut,buttonmask,hh)
%
% [answer] = EEGResponseBoxButton(TimeOut,buttonmask_,hh)
%
%
%  Copyright (C) 2017 Manfred Mauermann
%

answer     = 0;

% data   = guidata(at_hd);
% smpcfg = data.soundcfg;


% obtain actual rec-bufferlength and recpos
[success, recbuf, pos_old] = soundmexpro('recgetdata');
if success ~= 1,  soundmexpro(exit); error(['error calling ''recgetdata''' error_loc(dbstack)]);   end
Lrb = length(recbuf);

tic;
pause(0.001);
while ~answer & toc < TimeOut
    
    % check recording status 
    [success, value] = soundmexpro('recording');
    if 1 ~= success, soundmexpro('exit'); error('Fatal error in loop. Aborting...'); end
    
    % obtain actual rec-bufferlength and recpos
    [success, recbuf, pos_new] = soundmexpro('recgetdata');
    if success ~= 1,  soundmexpro(exit); error(['error calling ''recgetdata''' error_loc(dbstack)]);   end

    if pos_new+Lrb < 1, continue; end       
    
    diffpos = pos_new - pos_old;
    overlap = Lrb - diffpos;
    
    if ~diffpos, continue; end   
    
    if diffpos >= Lrb 
        warning('some response samples were missed!')
        overlap = 1;
    end

    pos_old = pos_new;
    %=================================
    relevantrecbuf = recbuf(max(overlap,1):end);% new samples
%     
%     figure(1)
%     plot(relevantrecbuf)
    %=================================
    
    DiffIDx   = find(diff(relevantrecbuf),1)+1;
    if ~isempty(DiffIDx)        
        amplitude      = (relevantrecbuf(DiffIDx)-relevantrecbuf(DiffIDx-1));
        press_release  = sign(amplitude);%- from buttonpressed to button released, + vice verca
        button         = round(log2(abs(amplitude).*2.^16));
        answer         = buttonmask(button);
        answer         = answer(press_release > 0);
        %answer_press   = answer(press_release > 0);
        %answer_release = answer(press_release < 0);
        if isempty(answer), answer = 0;continue; end
        if length(answer) > 1
            answer = answer(1);%erste positive (press) Antwort
            warning('multiple buttons has been pressed at once - first pressed keynumber will be used here')
        end
    end
    
    %check for controlbuttons via GUI like pause etc. data.answer will be
    %set anytime by the GUI-Callback functions
    drawnow;
    data = guidata(hh);
    data.answer
    if data.answer < 0
        answer = data.answer;
        break;
    end


end

