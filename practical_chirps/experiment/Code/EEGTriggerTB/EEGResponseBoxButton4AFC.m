%------------------------------------------------------------------------------
% EEG Trigger-Toolbox Version 0.2
%
% Author(s): Manfred Mauermann
%
% Copyright (c) 2017-2018, Manfred Mauermann. 
% All rights reserved.
%
% This work is licensed under the 
% Creative Commons Attribution-NonCommercial-NoDerivs 4.0 International License (CC BY-NC-ND 4.0). 
% To view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative
% Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
%------------------------------------------------------------------------------
%
% last modified 17.08.2018 Manfred Mauermann
%

function [answer,y] = EEGResponseBoxButton4AFC
%
% only to be used with an modified version 14 of the 
% AFC for Mathwork’s MATLAB (Version 1.40.0) by Stephan Ewert.
%
% The [answer,y] = EEGResponseBoxButton4AFC function which is waiting 
% for a Responsebox input and returning the "answer". 
% The function will stop and return  the answer if the answer isgiven via 
% the EEG-Responsebox or will stop only if the answer is given with 
% another device (keyboard, mouse, Touchscreen) via AFC. 
%


% flags used in AFC typically as variable y to describe the input status of the GUI 
%  0 = not ready, 
%  1 = accept commands,
%  2 = accept only end and restart command,
% -1 = proceed in afc_work if terminate and end flags are not set
% -2 = accept only end
%  4 = any keyboard and window button (waitfor in afc_main)

global eeg
global work

h = findobj('Tag','afc_win');
y = get(h,'UserData');

if y == 0 % if the AFC-GUI is not ready for input, 
          %the EEG-Input-Box should not be ready as well
    return;
end

TimeOut = inf;

set(h,'UserData',1)
set(h,'currentcharacter',char(0));

% obtain first rec-bufferlength and recpos
[success, recbuf, pos_old] = soundmexpro('recgetdata');
if success ~= 1,  soundmexpro(exit); error(['error calling ''recgetdata''' error_loc(dbstack)]);   end
Lrb = length(recbuf);


tic;
pause(0.001);

%============ wait/check for keypress on the EEG-response box ================
set(h,'UserData',1);
answer = 0;

% in the while loop soundmex is listening to the SPDIF-Input channel
% if a  Response-Box button specific bit is set. Note these bits are independend
% from the trigger bits. The respevtive answer will be given back.
% The loop will run until either a answer is given via Response-Box an is
% registered within this loop or an answer is given via the GUI-Callback
% function "afc_pressfnc.m", in this case "afc_pressfnc.m" will set the
% UserData of 'afc_win' and y will be ~=1 ot the work.abortall flag is set.
% or an timeout occurs (actuallly this is set to inf).
while ~answer & toc < TimeOut & ~work.abortall & y==1
    y = get(h,'UserData');
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
    %=================================
    DiffIDxOn = []; DiffIDxOff = [];

    DiffIDxOn   = find(diff(relevantrecbuf),1)+1;
    DiffIDxOff  = find(diff(relevantrecbuf),1)-1;
    
    try
    if ~isempty(DiffIDxOn)
        onset  = abs(relevantrecbuf(DiffIDxOn))./2;
    end
    if ~isempty(DiffIDxOff)
        offset = abs(relevantrecbuf(DiffIDxOff))./2;
    end
    catch
       warning('urgh DiffIDxOn: %3.1f   DiffIDxOff: %3.1f',DiffIDxOn,DiffIDxOff);
       keyboard 
    end
    
    if ~isempty(DiffIDxOn) |  ~isempty(DiffIDxOff)
        resp = max([onset,offset]).*2.^16;% these values are independent from the EEG system
        if length(resp) > 1 %maybe meanwhile obsolete
            warning('multiple buttons has been pressed at once - lowest pressed keynumber will be used here');
        end
        answer  = eeg.buttonmask(log2(resp(1))+1);
        %direction = sign(ceil(onset)-ceil(offset));           
        %if direction == -1
        %    answer = 0; % only trigger onset is relevant here
        %    continue;
        %end
    end
    drawnow;
end



