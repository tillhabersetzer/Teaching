% settings

% ensure that we don't mix up settings
clear settings

% Define paths
%-------------
settings.path2fieldtrip = fullfile('C:','Users','tillhabersetzer','Nextcloud','Synchronisation','Projekte','Toolboxen','fieldtrip-20221009');
% settings.path2fieldtrip = fullfile('C:','Users','tillhabersetzer','Nextcloud','Synchronisation','Projekte','Toolboxen','fieldtrip-20220917');
settings.path2project   = fullfile('M:','Blockpraktikum2022','analysis');

% settings.path2openmeeg = fullfile('C:','Users','tillhabersetzer','Nextcloud','Synchronisation','Projekte','Toolboxen','OpenMEEG');
% settings.path2iso2mesh = fullfile('C:','Users','tillhabersetzer','Nextcloud','Synchronisation','Projekte','Toolboxen','iso2mesh');
% settings.path2iso2cpd2 = fullfile('C:','Users','tillhabersetzer','Nextcloud','Synchronisation','Projekte','Toolboxen','cpd2');
% settings.path2dicm2nii = fullfile('C:','Users','tillhabersetzer','Nextcloud','Synchronisation','Projekte','Toolboxen','github_repo');

settings.conditions = {'clicks','upchirps','downchirps'};
settings.maxfilter = true; % apply maxfilter in data analysis
settings.latency_correction = 5/1000; % 5ms
