function hfig = display_fixation_cross()
% to handle pos = get(0, 'MonitorPosition') right, all monitors must be
% connected before the start of matlab

    monitor = 1; % choose monitor for the fixation cross screen
    imshow('Fixationskreuz.png')
    hfig = gcf;
    pos = get(0, 'MonitorPosition');
    set(hfig, 'Position', pos(monitor,:),'Color','black','name','Fixationskreuz'); % set position of first monitor (left one)
    set(hfig,'WindowState','fullscreen','menubar','none','NumberTitle','off');

end