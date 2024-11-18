close all
clear 
clc 

% https://www.fieldtriptoolbox.org/workshop/natmeg/preprocessing/#mismatch-negativity
% Calculate scalp-current density with ft_scalpcurrentdensity????
% create eeg layout

% https://www.fieldtriptoolbox.org/tutorial/continuous/
% https://www.fieldtriptoolbox.org/tutorial/preprocessing_erp/

% avg reference better than nose reference

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

% acoustic latency correction for triggers
latency_correction = settings.latency_correction;

conditions = settings.conditions;

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

%% Compute N19P30 32 channel EEG CAP
%--------------------------------------------------------------------------
subidx  = 1; % has only been measured for sub-01
subject = ['sub-',num2str(subidx,'%02d')];

% layout
%--------------------------------------------------------------------------
datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-clicks.fif']);
elec     = ft_read_sens(datapath, 'senstype', 'eeg');

cfg         = [];
cfg.elec    = elec;
cfg.channel = 1:32; % 32 EEG channels
layout      = ft_prepare_layout(cfg);

% cfg        = [];
% cfg.layout = layout;
% ft_layoutplot(cfg)

% Preprocess data
%--------------------------------------------------------------------------

avgs     = cell(2,3);
N_trials = zeros(2,3);

for ridx = 1:2 % loop over reference

    for cidx = 1:3 % loop over conditions
    
        % Path for EEG data
        %------------------
        condition = conditions{cidx};
        datapath  = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);
    
        switch condition
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
    
        % Filter continuous data to avoid edge artifacts
        %-----------------------------------------------
        cfg              = [];
        cfg.dataset      = datapath;
        cfg.chantype     = 'EEG';
        cfg.channel      = 4:35; % 32 EEG channels
        cfg.continuous   = 'yes';
        cfg.bpfilter     = 'yes';
        cfg.bpfreq       = [16,120];
        cfg.dftfilter    = 'yes';       % enable notch filtering to eliminate power line noise
        cfg.dftfreq      = [50,100];    % set up the frequencies for notch filtering
        cfg.dftreplace   = 'zero';      % 'zero' is default'
        cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
        data             = ft_preprocessing(cfg); 
    
        % rereference
        % Meiser: nose-tip reference
        if ridx==2 % avg reference
            cfg            = [];
            cfg.reref      = 'yes';
            cfg.refmethod  = 'avg';
            cfg.channel    = 'all'; % this is the default
            cfg.refchannel = 'all'; % average is computed over the channels specified in cfg.refchannel
            data           = ft_preprocessing(cfg,data);
        end
    
        % Epoch data
        %-----------
        cfg     = [];
        cfg.trl = trl;  
        epoched = ft_redefinetrial(cfg,data); 
    
        % Demean trials
        %--------------
        cfg                = [];
        cfg.demean         = 'yes';
        cfg.baselinewindow = [-0.05 0];
        epoched            = ft_preprocessing(cfg,epoched);  
    
        % Semi automatic artifact rejection
        %----------------------------------
%         cfg        = [];
%         cfg.metric = 'zvalue';
%         epoched    = ft_rejectvisual(cfg,epoched);

        % Note down amount of preserved epochs
        %-------------------------------------
        N_trials(ridx,cidx) = length(epoched.trial);
    
        % Timelockanalysis
        %-----------------
        cfg = [];
        avg = ft_timelockanalysis(cfg, epoched);
         % Add latency correction 
        %-----------------------
        % ~5ms 
        % 1.8ms transmission delay of sound + max 3ms Recording of analog MEG signals
        % -> digital triggers arrive 5 ms to early in MEG recordings
        avg.time = avg.time - latency_correction; % in sec

        avgs{ridx,cidx} = avg;
        clear data epoched avg
    
    end % conditions

end % electrode references

% Save data
%--------------------------------------------------------------------------
reference = {'nose-tip reference','avg reference'};
% make folder for data
dir2save = fullfile(settings.path2project,'derivatives',subject,'sensorlevel');
if ~exist(dir2save,'dir')
    mkdir(dir2save)
end
save(fullfile(dir2save,[subject,'_erp-N19P30.mat']),'avgs','N_trials','reference','layout'); 

%% Compute N19P30 10 channel cEEGrid
%--------------------------------------------------------------------------

subjects = [1,2];

% Rereferencing/ channel combinations
%--------------------------------------------------------------------------
combinations    = cell(1,3);
% combinations left/ right ear
combinations{1} = nchoosek(1:5,2); % right ear
combinations{2} = nchoosek(6:10,2); % left ear
% combinations between ears
[m,n]           = meshgrid(1:5,6:10);
combinations{3} = [m(:),n(:)];

combinations_label    = cell(1,3);
combinations_label{1} = cell(size(combinations{1}));
combinations_label{2} = cell(size(combinations{2}));
combinations_label{3} = cell(size(combinations{3}));

% left / right ear
for c = 1:numel(combinations{1}) 
    combinations_label{1}{c} = ['EEG0',num2str(combinations{1}(c),'%02d')]; % right ear
    combinations_label{2}{c} = ['EEG0',num2str(combinations{2}(c),'%02d')]; % left ear
end
% combinations between ears
for c = 1:numel(combinations{3})
    combinations_label{3}{c} = ['EEG0',num2str(combinations{3}(c),'%02d')]; % right ear
end

C(1)        = size(combinations{1},1);
C(2)        = size(combinations{3},1);
avgs_left    = cell(C(1),3);
avgs_right   = cell(C(1),3);
avgs_between = cell(C(2),3);

% Computations over subjects
%--------------------------------------------------------------------------
for subidx=subjects % loop over subjects
    subject = ['sub-',num2str(subidx,'%02d')];
    
    % Preprocess EEG
    %----------------------------------------------------------------------
    data = cell(1,3);
    trl  = cell(1,3);
    
    for cidx = 1:3 % loop over conditions
        
        % Path for EEG data
        %------------------
        condition = conditions{cidx};
        datapath  = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);
    
        switch condition
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
        trl{cidx}               = cfg.trl;
    
        % Filter continuous data to avoid edge artifacts
        %-----------------------------------------------
        cfg              = [];
        cfg.dataset      = datapath;
        cfg.channel      = 'EEG';
        cfg.continuous   = 'yes';
        cfg.bpfilter     = 'yes';
        cfg.bpfreq       = [16,120];
        cfg.dftfilter    = 'yes';       % enable notch filtering to eliminate power line noise
        cfg.dftfreq      = [50,100];    % set up the frequencies for notch filtering
        cfg.dftreplace   = 'zero';      % 'zero' is default'
        cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
        data{cidx}       = ft_preprocessing(cfg); 
    end
    
    % Rename channels
    %----------------------------------------------------------------------
    % sub-01: 32 channel EEG + 10 channel cEEGrid
    % sub-02: 10 channel cEEGrid
    if strcmp(subject,'sub-01')
        for cidx = 1:3
            cfg         = [];
            cfg.channel = 33:42;
            data{cidx}  = ft_selectdata(cfg,data{cidx});
            % Rename channels
            for c = 1:10
                data{cidx}.label{c} = ['EEG0',num2str(c,'%02d')];
            end
        end
    end
    
    % layout
    %--------------------------------------------------------------------------
    cfg      = [];
    cfg.elec = data{1}.elec;
    if strcmp(subject,'sub-01')
        cfg.channel = 33:42;
    end
    if strcmp(subject,'sub-02')
        cfg.channel = 1:10;
    end
    layout      = ft_prepare_layout(cfg);

    % rename channels
    if strcmp(subject,'sub-01')
        for c = 1:10
            layout.label{c} = ['EEG0',num2str(c,'%02d')];
        end
    end
        
    cfg        = [];
    cfg.layout = layout;
    ft_layoutplot(cfg)
    
    % Referencing / Computation of bipolar channels
    %----------------------------------------------------------------------  
    for cidx = 1:3 % loop over conditions

        for eidx = 1:3 % loop over ears (left / right / between)

            switch eidx
                case {1,2}
                    C_num = C(1);
                case 3
                    C_num = C(2);
            end

            for c = 1:C_num % loop over bipolar channels

                % rereference
                %------------
                cfg            = [];
                cfg.channel    = combinations_label{eidx}(c,:);
                cfg.reref      = 'yes';
                cfg.refchannel = combinations_label{eidx}{c,2};
                data1          = ft_preprocessing(cfg,data{cidx});
        
                cfg         = [];
                cfg.channel = combinations_label{eidx}{c,1};
                data1       = ft_preprocessing(cfg, data1);
            
                % Epoch data
                %-----------
                cfg      = [];
                cfg.trl  = trl{cidx};  
                epoched1 = ft_redefinetrial(cfg,data1); 
                
                % Demean trials
                %--------------
                cfg                = [];
                cfg.demean         = 'yes';
                cfg.baselinewindow = [-0.05 0];
                epoched1           = ft_preprocessing(cfg,epoched1);  
                
                % Semi automatic artifact rejection
                %----------------------------------
            %     cfg        = [];
            %     cfg.metric = 'zvalue';
            %     epoched1   = ft_rejectvisual(cfg,epoched1);
            
                % Timelockanalysis
                %-----------------
                cfg = [];
                avg = ft_timelockanalysis(cfg, epoched1);

                % Add latency correction 
                %-----------------------
                % ~5ms 
                % 1.8ms transmission delay of sound + max 3ms Recording of analog MEG signals
                % -> digital triggers arrive 5 ms to early in MEG recordings
                avg.time = avg.time - latency_correction; % in sec
              
                if eidx==1
                    avgs_right{c,cidx} = avg;
                elseif eidx==2
                    avgs_left{c,cidx} = avg;
                elseif eidx==3
                    avgs_between{c,cidx} = avg;
                end 

                clear data1 epoched1 avg

            end

        end
    end
    
    
    % Save data
    %--------------------------------------------------------------------------
    % make folder for data
    dir2save = fullfile(settings.path2project,'derivatives',subject,'sensorlevel');
    if ~exist(dir2save,'dir')
        mkdir(dir2save)
    end
    save(fullfile(dir2save,[subject,'_erp-N19P30_cEEGrid.mat']),'avgs_left','avgs_right','avgs_between','combinations_label','layout'); 

    clear avgs_left avgs_right avgs_between layout

end

%% Missmatch negativity - cEEGrid
%--------------------------------------------------------------------------
subidx  = 2; % has only been measured for sub-01
subject = ['sub-',num2str(subidx,'%02d')];

% Rereferencing/ channel combinations
%--------------------------------------------------------------------------
combinations    = cell(1,3);
% combinations left/ right ear
combinations{1} = nchoosek(1:5,2); % right ear
combinations{2} = nchoosek(6:10,2); % left ear
% combinations between ears
[m,n]           = meshgrid(1:5,6:10);
combinations{3} = [m(:),n(:)];

combinations_label    = cell(1,3);
combinations_label{1} = cell(size(combinations{1}));
combinations_label{2} = cell(size(combinations{2}));
combinations_label{3} = cell(size(combinations{3}));

% left / right ear
for c = 1:numel(combinations{1}) 
    combinations_label{1}{c} = ['EEG0',num2str(combinations{1}(c),'%02d')]; % right ear
    combinations_label{2}{c} = ['EEG0',num2str(combinations{2}(c),'%02d')]; % left ear
end
% combinations between ears
for c = 1:numel(combinations{3})
    combinations_label{3}{c} = ['EEG0',num2str(combinations{3}(c),'%02d')]; % right ear
end

C(1)         = size(combinations{1},1);
C(2)         = size(combinations{3},1);
avgs_left    = cell(C(1),3);
avgs_right   = cell(C(1),3);
avgs_between = cell(C(2),3);

% Computations 
%--------------------------------------------------------------------------

% Path for EEG data
datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-mmn.fif']);

cfg              = [];
cfg.dataset      = datapath;
cfg.channel      = 'EEG'; 
cfg.continuous   = 'yes';
cfg.bpfilter     = 'yes';
cfg.bpfreq       = [1,30];
cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
data             = ft_preprocessing(cfg); 

% Define trials
%--------------
cfg                     = [];
cfg.dataset             = datapath;
cfg.trialfun            = 'ft_trialfun_general'; 
cfg.trialdef.eventtype  = 'STI101';
cfg.trialdef.eventvalue = [1,2]; % standard / deviant
cfg.trialdef.prestim    = 0.2;                  
cfg.trialdef.poststim   = 0.8;                 
cfg                     = ft_definetrial(cfg);
trl                     = cfg.trl;

% layout
%--------------------------------------------------------------------------
cfg         = [];
cfg.elec    = data.elec;
cfg.channel = 1:10;
layout      = ft_prepare_layout(cfg);

cfg        = [];
cfg.layout = layout;
ft_layoutplot(cfg)

% Referencing / Computation of bipolar channels
%----------------------------------------------------------------------  
for eidx = 1:3 % loop over ears (left / right / between)

    switch eidx
        case {1,2}
            C_num = C(1);
        case 3
            C_num = C(2);
    end

    for c = 1:C_num % loop over bipolar channels

        % rereference
        %------------
        cfg            = [];
        cfg.channel    = combinations_label{eidx}(c,:);
        cfg.reref      = 'yes';
        cfg.refchannel = combinations_label{eidx}{c,2};
        data1          = ft_preprocessing(cfg,data);

        cfg         = [];
        cfg.channel = combinations_label{eidx}{c,1};
        data1       = ft_preprocessing(cfg, data1);
    
        % Epoch data
        %-----------
        cfg      = [];
        cfg.trl  = trl;  
        epoched1 = ft_redefinetrial(cfg,data1); 
        
        % Demean trials
        %--------------
        cfg                = [];
        cfg.demean         = 'yes';
        cfg.baselinewindow = [-0.2 0];
        epoched1           = ft_preprocessing(cfg,epoched1);  
        
        % Semi automatic artifact rejection
        %----------------------------------
    %     cfg        = [];
    %     cfg.metric = 'zvalue';
    %     epoched1   = ft_rejectvisual(cfg,epoched1);

        % Timelockanalysis
        %-----------------
        % standard tones
        cfg        = [];
        cfg.trials = find(epoched1.trialinfo(:,1) == 1);
        avg_sta    = ft_timelockanalysis(cfg, epoched1);
        
        % deviant tones
        cfg        = [];
        cfg.trials = find(epoched1.trialinfo(:,1) == 2);
        avg_dev    = ft_timelockanalysis(cfg, epoched1);

        % differences
        cfg           = [];
        cfg.operation = 'subtract';
        cfg.parameter = 'avg';
        cfg.operation = 'x1 - x2'; % see https://www.fieldtriptoolbox.org/workshop/natmeg/dipolefitting/
        avg_diff      = ft_math(cfg, avg_dev, avg_sta);

        % Add latency correction 
        %-----------------------
        % ~5ms 
        % 1.8ms transmission delay of sound + max 3ms Recording of analog MEG signals
        % -> digital triggers arrive 5 ms to early in MEG recordings
        avg_sta.time  = avg_sta.time - latency_correction; % in sec
        avg_dev.time  = avg_dev.time - latency_correction; % in sec
        avg_diff.time = avg_diff.time - latency_correction; % in sec
      
        if eidx==1
            avgs_right{c,1} = avg_sta;
            avgs_right{c,2} = avg_dev;
            avgs_right{c,3} = avg_diff;
        elseif eidx==2
            avgs_left{c,1} = avg_sta;
            avgs_left{c,2} = avg_dev;
            avgs_left{c,3} = avg_diff;
        elseif eidx==3
            avgs_between{c,1} = avg_sta;
            avgs_between{c,2} = avg_dev;
            avgs_between{c,3} = avg_diff;
        end 

        clear data1 epoched1 avg_sta avg_dev avg_diff

    end
end

order = {'standard','deviant','difference'};

% Save data
%--------------------------------------------------------------------------
% make folder for data
dir2save = fullfile(settings.path2project,'derivatives',subject,'sensorlevel');
if ~exist(dir2save,'dir')
mkdir(dir2save)
end
save(fullfile(dir2save,[subject,'_erp-mmn_cEEGrid.mat']),'avgs_left','avgs_right','avgs_between','combinations_label','layout','order'); 

clear avgs_left avgs_right avgs_between layout

%% Clean-Up
rmpath('..')

