% settings

% ensure that we don't mix up settings
clear settings

% Define paths
%-------------
settings.path2fieldtrip = fullfile('C:\Users\tillhabersetzer\Nextcloud\Synchronisation\Projekte\Toolboxen\fieldtrip-20240620');
settings.path2project   = fullfile('E:\practical_chirps_v2\practical_meg');

settings.conditions         = {'clicks','upchirps','downchirps'};
settings.maxfilter          = false; % apply maxfilter in data analysis
settings.latency_correction = 5/1000; % 5ms
