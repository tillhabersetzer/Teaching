%--------------------------------------------------------------------------
% Plot sensorlevel results
%
% Check out: https://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging/
%--------------------------------------------------------------------------

close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','subjectdata'))
eval('main_settings')

addpath(fullfile('..','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
subject = 'sub-04';

% Choose which type of evoked field (erf) to fit
% erf_type = 'N19mP30m'; % type for practical
erf_type = 'N100m'; 

% Choose channel type
% chantype = 'meg'; 
chantype = 'eeg'; 

% order in avg cell
stimtype = {'clicks','upchirps','downchirps'};
%--------------------------------------------------------------------------


%% Load data
%-----------

% Switch between erf-types
%-------------------------
switch erf_type
    case 'N19mP30m'
        data    = importdata(fullfile(settings.path2bids,'derivatives',subject,'sensorlevel',[subject,'_erf-N19mP30m_',chantype,'.mat']));    
        timewin = [0.03 0.05];
    case 'N100m'
        data = importdata(fullfile(settings.path2bids,'derivatives',subject,'sensorlevel',[subject,'_erf-N100m_',chantype,'.mat']));
        timewin = [0.095,0.125];
end 

avg      = data.avg;
N_trials = data.N_trials;

switch chantype
    case 'meg'
        layout = 'neuromag306mag.lay';
%         layout = 'neuromag306planar.lay';
        chan2plot = 'MEG1411';
    case 'eeg'
        layout = importdata(fullfile(settings.path2bids,'derivatives',subject,'layout',[subject,'_layout-eeg_orig.mat']));   
        chan2plot = 'EEG023';
end
%% Plot of event related fields for all sensors arranged topographically 
%--------------------------------------------------------------------------
% plot all 3 conditions together

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
% cfg.layout     = 'neuromag306planar.lay';
cfg.linecolor  = 'rbk';
ft_multiplotER(cfg, avg{:});
legend({'clicks';'upchirps';'downchirps'});
title([chantype,': ',subject])

%% Plot one specific channel
%--------------------------------------------------------------------------
% plot all 3 conditions together

cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg{:});
legend({'clicks';'upchirps';'downchirps'});
title([chantype,': ',subject])

%% Plot the topographic distribution of the data
%--------------------------------------------------------------------------
% plot over conditions

for cidx=1:3
    cfg          = [];
    cfg.xlim     = timewin; 
    cfg.colorbar = 'yes';
    cfg.layout   = layout;
%     cfg.layout   = 'neuromag306planar.lay';
    ft_topoplotER(cfg, avg{cidx});
    title([chantype,': ', subject,' ',stimtype{cidx}])
end


%% Combine planar gradients of averaged data
%--------------------------------------------------------------------------
if strcmp(chantype,'meg')
    avgplanarComb = cell(1,3);
    
    for cidx=1:3
        cfg                 = [];
        avgplanarComb{cidx} = ft_combineplanar(cfg, avg{cidx});
    end
    
    
    % here: all conditions for one subject (s is fixed, loop over c)
    for cidx=1:3
        cfg        = [];
        cfg.xlim   = timewin; % average over 30-50 ms
        cfg.layout = 'neuromag306cmb.lay';
        ft_topoplotER(cfg, avgplanarComb{cidx});
        title([subject,' ',stimtype{cidx}])
    end
end
