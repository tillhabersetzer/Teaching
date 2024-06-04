close all
clear all
clc

% load channel labels for channel selection
label = importdata('meg_label.mat');

% Step by step description for the Neuromag acquisition computer
%--------------------------------------------------------------------------
% 0.) Download and unzip fieldtrip

% 1.) run in a Linux terminal:
% [meg@sinuhe ~]$ . fieldtrip_rt.sh

% 2.) Start Acquisition. You should see some information being printed 
%     in the terminal that you used to start neuromag2ftx

% Step by step description for the visualization computer
%--------------------------------------------------------------------------

% 0.) Download and unzip fieldtrip

% 1.) Start MATLAB and make sure the online_meg directory is added to path:
% addpath 'C:\Program Files\MATLAB\fieldtrip-20191021'
% addpath(['D:' filesep 'Experiments' filesep 'asap' filesep 'fieldtrip-master']) 
% ft_defaults

% 2.) Type in MATLAB command window
cfg         = [];
cfg.dataset = 'buffer://134.106.116.42:1972';  % NEED TO CHANGE TO CORRECT IP HERE!
%ft_realtime_headlocalizer(cfg);               % does not work

cfg.blocksize = 1;                             % (default = 1 second)
% choose desired channels:
desired = {'MF-M1'}; % frontral magnetometer
%desired = {'MF-MG'}; % frontral gradiometer
%desired = {'M1'};    % only Magnetometer

my_selection = channel_selection(desired,label);
cfg.channel  = my_selection;
ft_realtime_signalviewer(cfg)

%--------------------------------------------------------------------------
% press Ctrl-C to stop the realtime function (in Linux Terminal + MATLAB)
%--------------------------------------------------------------------------




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


