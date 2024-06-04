%--------------------------------------------------------------------------
% Computation of event related fields over subjects (1,2,3,4) and 
% conditions {'clicks','upchirps','downchirps'}
%
% Check out: https://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging/
%
% Up to now:
% - ignore artifact channels (EOG, ECG), no ICA, no SSP applied
% - use simple z-value threshold for rejection of bad epochs, might be too 
%   easy and naiv
% - no headposition transformation and movement correction applied yet,
%   only simple tSSS method. The first point in particular may lead to 
%   questionable results when comparing different conditions because the
%   conditions were recorded in separate runs and hence have different head
%   positions.
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
% subjects = [1,2,3,4];
subjects = 4;

stimtype = {'clicks','upchirps','downchirps'};

% chantype = {'meg','eeg'}; 
% chantype = {'meg'};
chantype = {'eeg'};
% chantype = {'megmag'};
%--------------------------------------------------------------------------

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

Nsub  = length(subjects);
Nsens = length(chantype);

%% Preprocess data
%--------------------------------------------------------------------------

parfor subidx=subjects % loop over subjects

    data              = cell(1,3);
    data_preprocessed = cell(1,3);
    avg               = cell(1,3);

    for sidx=1:Nsens % loop over channel types

        N_trials = zeros(1,3);
    
        for cidx=1:3 % loop over conditions
    
            % Path for MEG data
            %------------------------------------------------------------------
            % Path for MEG data (no headposition transformation applied)
            subject  = ['sub-',num2str(subidx,'%02d')];
            suffix   = '-raw_tsss';
            datapath = fullfile(settings.path2bids,'derivatives',subject,'tsss',[subject,'_task-',stimtype{cidx},'_meg',suffix,'.fif']);
        
            switch stimtype{cidx}
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
            cfg.trialfun            = 'ft_trialfun_general'; % this is the default
            cfg.trialdef.eventtype  = 'STI101';
            cfg.trialdef.eventvalue = eventvalue; 
            cfg.trialdef.prestim    = 0.05;                  % in seconds
            cfg.trialdef.poststim   = 0.35;                  % in seconds
            cfg                     = ft_definetrial(cfg);
            trl                     = cfg.trl;

            switch chantype{sidx}
                case 'meg'
                    
                    % Filter continuous data to avoid edge artifacts
                    %-----------------------------------------------
                    cfg              = [];
                    cfg.dataset      = datapath;
                    cfg.channel      = 'meg'; 
                    cfg.continuous   = 'yes';
                    cfg.bpfilter     = 'yes';
                    cfg.bpfreq       = [4,30];
                    cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
                    data             = ft_preprocessing(cfg); 

                    % Remove bad trials
                    %------------------
                    % Jump
                    cfg                        = [];
                    cfg.trl                    = trl;
                    cfg.dataset                = datapath;
                    cfg.artfctdef.jump.channel = 'meg'; 
                    [~, artifacts_jump]        = ft_artifact_jump(cfg);
        
                    % Clip
                    cfg                        = [];
                    cfg.trl                    = trl;
                    cfg.dataset                = datapath;
                    cfg.artfctdef.clip.channel = 'meg';
                    [~, artifacts_clip]        = ft_artifact_clip(cfg);
 
                    % Muscle
                    cfg                          = [];
                    cfg.trl                      = trl;
                    cfg.dataset                  = datapath;
                    cfg.artfctdef.muscle.channel = 'meg';
            %         cfg.artfctdef.muscle.interactive = 'yes';
                    cfg.artfctdef.muscle.cutoff  = 8; % default 4
                    [~, artifacts_muscle]        = ft_artifact_muscle(cfg);
        
                    % Update trl-matrix
                    cfg                           = [];
                    cfg.trl                       = trl;
                    cfg.dataset                   = datapath;
                    cfg.artfctdef.jump.artifact   = artifacts_jump;
                    cfg.artfctdef.clip.artifact   = artifacts_clip;
                    cfg.artfctdef.muscle.artifact = artifacts_muscle;
                    cfg                           = ft_rejectartifact(cfg);
                    
                    trl_new = cfg.trl;

                case 'eeg'

                    % Filter continuous data to avoid edge artifacts
                    %-----------------------------------------------
                    cfg              = [];
                    cfg.reref        = 'yes';
                    cfg.refchannel   = 'all';
                    cfg.dataset      = datapath;
                    cfg.channel      = 'eeg'; 
                    cfg.continuous   = 'yes';
                    cfg.bpfilter     = 'yes';
                    cfg.bpfreq       = [4,30];
                    cfg.coilaccuracy = 0;   
                    cfg.detrend      = 'yes';
                    data             = ft_preprocessing(cfg); 

                    % Remove bad trials
                    %------------------  
                    % Clip
                    cfg                        = [];
                    cfg.trl                    = trl;
                    cfg.dataset                = datapath;
                    cfg.artfctdef.clip.channel = 'eeg';
                    [~, artifacts_clip]        = ft_artifact_clip(cfg);
 
                    % Muscle
                    cfg                          = [];
                    cfg.trl                      = trl;
                    cfg.dataset                  = datapath;
                    cfg.artfctdef.muscle.channel = 'eeg';
%                     cfg.artfctdef.muscle.interactive = 'yes';
                    cfg.artfctdef.muscle.cutoff  = 8; % default 4
                    [~, artifacts_muscle]        = ft_artifact_muscle(cfg);
        
                    % Update trl-matrix
                    cfg                           = [];
                    cfg.trl                       = trl;
                    cfg.dataset                   = datapath;
                    cfg.artfctdef.clip.artifact   = artifacts_clip;
                    cfg.artfctdef.muscle.artifact = artifacts_muscle;
                    cfg                           = ft_rejectartifact(cfg);
                    
                    trl_new = cfg.trl;

            end
                
                % Epoch data
                %-----------
                cfg     = [];
                cfg.trl = trl_new;  
                data    = ft_redefinetrial(cfg,data); 

                % Demean trials
                %--------------
                cfg                = [];
                cfg.demean         = 'yes';
                cfg.baselinewindow = [-0.05 0];
                data               = ft_preprocessing(cfg,data);  

                % Note down amount of preserved epochs
                %-------------------------------------
                N_trials(cidx) = length(data.trial);     

                % Rename data for pooling
                %------------------------
                data_preprocessed{cidx} = data;

                %% Timelockanalysis
                cfg       = [];
                avg{cidx} = ft_timelockanalysis(cfg, data);

        end % loop conditions

        % Compute pooled average (over clicks and downchirps)
        %----------------------------------------------------

        cfg              = [];
        data_pooled      = ft_appenddata(cfg, data_preprocessed{1}, data_preprocessed{3});
        data_pooled.grad = data_preprocessed{1}.grad;
        if isfield(data_preprocessed,'elec')
            data_pooled.elec = data_preprocessed{1}.elec;
        end
  
        cfg        = [];
        avg_pooled = ft_timelockanalysis(cfg, data_pooled);

        % Save data
        %----------
        % make folder for data
        dir2save = fullfile(settings.path2bids,'derivatives',subject,'sensorlevel');
        if ~exist(dir2save, 'dir')
            mkdir(dir2save)
        end
        parsave(fullfile(dir2save,[subject,'_erf-N100m_',chantype{sidx},'.mat']), ...
            avg,avg_pooled,N_trials); 

%         clear avg avg_pooled N_trials data_preprocessed
   
    end % loop sensor types
 
end % loop subjects

%% functions
%--------------------------------------------------------------------------
function parsave(fname,avg,avg_pooled,N_trials)
  save(fname,'avg','avg_pooled','N_trials')
end
