%--------------------------------------------------------------------------
% Computation of event related fields over subjects (1,2,...) and 
% conditions {'clicks','upchirps','downchirps'}
%
% Check out: https://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging/
% for additional information
%
% Up to now:
% - Artifact channels (HEOG, VEOG, ECG) are not used in analysis. 
%   Application of ICA (or SSP - not supported in fieldtrip yet) is 
%   possible but ignored due to increased complexity in pipeline. 
%   Furthermore, the benefit of these artifact suppression methods seems 
%   to be fairly small.
% - Instead contaminated epochs are rejected with a semi-automatic
%   procedure based on calculated z-values. Epochs above an individually
%   determined threshold are discarded. 
% - If maxfiltered data is analyzed, the MEG recordings of different 
%   conditions (clicks, upchirps, downchirps) are optionally transformed
%   to a common headposition between conditions so that a common forward
%   model could be used during the source space analysis (dipole fitting).
%   No application of a head position transformation may lead to 
%   questionable results when comparing different conditions because the
%   conditions were recorded in separate runs and therefore have different 
%   head positions.
%--------------------------------------------------------------------------

close all
clear 
clc 

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

% check whether maxfiltered data should be analyzed
maxfilter = settings.maxfilter;
% acoustic latency correction for triggers
latency_correction = settings.latency_correction;

%% Script settings
%--------------------------------------------------------------------------
subjects  = [1,2];
%--------------------------------------------------------------------------

% conditions
conditions = settings.conditions;

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

%% Preprocess data
%--------------------------------------------------------------------------

for subidx=subjects % loop over subjects

    subject           = ['sub-',num2str(subidx,'%02d')];
    data_preprocessed = cell(1,3);
    avg               = cell(1,3);
    N_trials          = zeros(1,3);

    for cidx=1:3 % loop over conditions

        % Path for MEG data
        %------------------------------------------------------------------
        condition = conditions{cidx};
   
        if maxfilter
            datapath = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-',condition,'-raw_tsss.fif']);
        else
            datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);
        end

        switch conditions{cidx}
            case 'clicks'
                eventvalue = TrigID.click;
            case 'upchirps'
                eventvalue = TrigID.upchirp;
            case 'downchirps'
                eventvalue = TrigID.downchirp;
        end

        % Define trials
        %--------------
        cfg                     = [];
        cfg.dataset             = datapath;
        cfg.trialfun            = 'ft_trialfun_general'; 
        cfg.trialdef.eventtype  = 'STI101';
        cfg.trialdef.eventvalue = eventvalue; 
        cfg.trialdef.prestim    = 0.05;                  % in seconds
        cfg.trialdef.poststim   = 0.35;                  % in seconds
        cfg                     = ft_definetrial(cfg);
        trl                     = cfg.trl;
                
        % Filter continuous data to avoid edge artifacts
        %-----------------------------------------------
        cfg              = [];
        cfg.dataset      = datapath;
        cfg.channel      = 'meg'; 
        cfg.continuous   = 'yes';
        cfg.bpfilter     = 'yes';
        cfg.bpfreq       = [16,120];
        cfg.dftfilter    = 'yes';       % enable notch filtering to eliminate power line noise
        cfg.dftfreq      = [50,100];    % set up the frequencies for notch filtering
        cfg.dftreplace   = 'zero';      % 'zero' is default'
        cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
        data             = ft_preprocessing(cfg); 

        % Epoch data
        %-----------
        cfg     = [];
        cfg.trl = trl;  
        epoched = ft_redefinetrial(cfg,data); 
    
        % Semi automatic artifact rejection
        %----------------------------------
        % separately for magnetometers
        cfg             = [];
        cfg.metric      = 'zvalue';
        cfg.channel     = 'megmag';
        cfg.keepchannel = 'yes';  % This keeps those channels that are not displayed in the data
        epoched         = ft_rejectvisual(cfg,epoched);
        
        % separately for gradiometers
        cfg.channel     = 'megplanar';
        epoched         = ft_rejectvisual(cfg,epoched);

        % Demean trials
        %--------------
        cfg                = [];
        cfg.demean         = 'yes';
        cfg.baselinewindow = [-0.05 0];
        epoched            = ft_preprocessing(cfg,epoched);  

        data_preprocessed{cidx} = epoched;

        % Note down amount of preserved epochs
        %-------------------------------------
        N_trials(cidx) = length(epoched.trial);

        % Timelockanalysis
        %-----------------
        cfg       = [];
        avg{cidx} = ft_timelockanalysis(cfg, epoched);

        clear data epoched

    end % loop over conditions

    % Compute pooled average (over clicks and downchirps)
    %----------------------------------------------------
    idx              = find(contains(conditions,{'clicks','downchirps'}));

    cfg              = [];
    data_pooled      = ft_appenddata(cfg, data_preprocessed{idx(1)}, data_preprocessed{idx(2)});
    data_pooled.grad = data_preprocessed{1}.grad;

    cfg        = [];
    avg_pooled = ft_timelockanalysis(cfg, data_pooled);

    %% Add latency correction 
    %--------------------------------------------------------------------------
    % ~5ms 
    % 1.8ms transmission delay of sound + max 3ms Recording of analog MEG signals
    % -> digital triggers arrive 5 ms to early in MEG recordings

    for cidx=1:3
        avg{cidx}.time = avg{cidx}.time - latency_correction; % in sec
    end
    avg_pooled.time = avg_pooled.time - latency_correction;

    % Save data
    %----------
    % make folder for data
    dir2save = fullfile(settings.path2project,'derivatives',subject,'sensorlevel');
    if ~exist(dir2save,'dir')
        mkdir(dir2save)
    end
    save(fullfile(dir2save,[subject,'_erf-N19mP30m.mat']),'avg','avg_pooled','N_trials'); 
 
end % loop subjects

%% Clean-up
rmpath('..')
