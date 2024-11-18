%--------------------------------------------------------------------------
% Till Habersetzer, 12.08.2024
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Step by step description 
%--------------------------------------------------------------------------
% https://www.fieldtriptoolbox.org/development/realtime/neuromag/
%
% 0.) Setup fieldtrip realtime streaming
%
% 1.) Run in a Linux terminal:
%     .\neurmag2ft
%     folder:
%     fieldtrip/realtime/src/acquisition/neuromag/bin/x86_64-pc-linux-gnu/
%
% 2.) Start Acquisition. You should see some information being printed 
%     in the terminal that you used to start neuromag2ft
% 
% 3.) Run this script and fetch data of stream.
%
%--------------------------------------------------------------------------
% press Ctrl-C to stop the realtime function (in Linux Terminal + MATLAB)
%--------------------------------------------------------------------------

close all
clear all
clc

% Add fieldtrip path
%-------------------
% Already included in startup.m
% otherwise:
% addpath('path2fieldtrip')
% ft_defaults

% Load channel labels for channel selection
label = importdata('meg_label.mat');

% Execute realtime streaming
%---------------------------
desired = {'MF-M1'}; % frontral magnetometer
% desired = {'MF-MG'}; % frontral gradiometer
% desired = {'M1'};    % only Magnetometer
my_selection = channel_selection(desired,label);

cfg           = [];
cfg.blocksize = 1;       
cfg.channel   = 'all';
% cfg.channel   = my_selection;
cfg.dataset   = 'buffer://10.2.164.47:1972'; % NEED TO CHANGE TO CORRECT IP HERE!
ft_realtime_signalviewer(cfg)

% more elaborate stuff
%--------------------------------------------------------------------------

% realtime_average
%-----------------
% cfg = [];
% cfg.dataset  = 'buffer://localhost:1972';
% desired      = {'MF-M1'};                        % frontral magnetometer
% my_selection = channel_selection(desired,label);
% cfg.channel  = my_selection;
% cfg.trialfun = 'ft_trialfun_general';
% cfg.trialdef.eventtype  = 'STI101';
% cfg.trialdef.eventvalue = 2; 
% cfg.trialdef.prestim    = 0.05;                  % in seconds
% cfg.trialdef.poststim   = 0.35;                  % in seconds
% ft_realtime_average(cfg)

% realtime_powerestimate
%----------------------
% cfg           = [];
% cfg.dataset   = 'buffer://localhost:1972';       % where to read the data from
% desired       = {'MF-M1'};                       % frontral magnetometer
% my_selection  = channel_selection(desired,label);
% cfg.channel   = my_selection;
% cfg.blocksize = 1;                               % seconds
% cfg.foilim    = [0 200];                         % frequency-of-interest limits, Hz
% ft_realtime_powerestimate(cfg)


