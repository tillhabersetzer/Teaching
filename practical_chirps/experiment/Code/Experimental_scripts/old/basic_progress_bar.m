function basic_progress_bar(f,duration,SampFreq,update)
% f must be f = uifigure;
% duration in samples and 
% update defines length of pause in seconds for update

%f = uifigure;
d = uiprogressdlg(f,...
                  'Title','Fortschritt der Geschichte',...
                  'Icon','info',...
                  'ShowPercentage','on',...
                  'Cancelable','on');
samples    = 0;
frac       = samples/duration;
addsamples = update*SampFreq;


    while (samples < duration)
        
        if d.CancelRequested
             break
        end
        
        d.Value = frac;
        pause(update);
        samples = samples + addsamples;   
        frac    = samples/duration; 
                
    end
   
d.Value = 1;
%close(f)
close(d)
end


            