close all
clear all


monitor = 2; % choose monitor for the fixation cross screen
imshow('Fixationskreuz.png')
hfig = gcf;
pos = get(0, 'MonitorPosition');
set(hfig, 'Position', pos(monitor,:),'Color','black','name','Fixationskreuz'); % set position of first monitor (left one)
set(hfig,'WindowState','fullscreen','menubar','none','NumberTitle','off');

% get the figure and axes handles
 % hFig = gcf;
 % hAx  = gca;
 % set the figure to full screen
 % set(hFig,'units','normalized','outerposition',[0 0 1 1]);
 % set the axes to full screen
 % set(hAx,'Unit','normalized','Position',[0 0 1 1]);
 % hide the toolbar
 %set(hFig,'menubar','none')
 % to hide the title
 %set(hFig,'NumberTitle','off');
 
%set(hFig,'WindowState','fullscreen','Color','black')

