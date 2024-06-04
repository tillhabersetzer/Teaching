close all
clear all
clc
% 1.) Start MATLAB and make sure the online_meg directory is added to path:
%addpath 'C:\Program Files\MATLAB\fieldtrip-20191021'
addpath(['D:' filesep 'Experiments' filesep 'asap' filesep 'fieldtrip-master']) 
ft_defaults

% 2.) Type in MATLAB command window
cfg         = [];
cfg.dataset = 'buffer://134.106.116.42:1972';  % NEED TO CHANGE TO CORRECT IP HERE!
ft_realtime_headlocalizer(cfg);