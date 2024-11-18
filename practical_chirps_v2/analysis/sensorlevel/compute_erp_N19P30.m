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

%% Script settings
%--------------------------------------------------------------------------
subidx  = 1;
subject = ['sub-',num2str(subidx,'%02d')];

conditions = settings.conditions;
%--------------------------------------------------------------------------

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

%% Create EEG layout
datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-clicks.fif']);
elec     = ft_read_sens(datapath, 'senstype', 'eeg');

cfg         = [];
cfg.elec    = elec;
cfg.channel = 1:32; % 32 EEG channels
layout      = ft_prepare_layout(cfg);

cfg        = [];
cfg.layout = layout;
ft_layoutplot(cfg)

%% Preprocess data
%--------------------------------------------------------------------------

avg      = cell(2,3);
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
        cfg        = [];
        cfg.metric = 'zvalue';
        epoched    = ft_rejectvisual(cfg,epoched);
    
        % Timelockanalysis
        %-----------------
        cfg            = [];
        avg{ridx,cidx} = ft_timelockanalysis(cfg, epoched);
    
        % Note down amount of preserved epochs
        %-------------------------------------
        N_trials(ridx,cidx) = length(avg{ridx,cidx}.cfg.previous.trials);
    
        clear data epoched
    end % conditions

end % electrode references

%% Visualize data

% topoplot
%---------
cfg            = [];
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, avg{1,1}, avg{1,2}, avg{1,3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('nose reference')

cfg            = [];
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, avg{2,1}, avg{2,2}, avg{2,3});
legend({conditions{1};conditions{2};conditions{3}});
subtitle('avg reference')

% single channel
%---------------
chan2plot = 'EEG005';
cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg{:});
legend({conditions{1};conditions{2};conditions{3}});
title(subject)

%% Bipolar channels ("cEEGrid")

% Path for EEG data
%------------------
condition = 'upchirps';
datapath  = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);

% layout
%--------------------------------------------------------------------------
elec     = ft_read_sens(datapath, 'senstype', 'eeg');

cfg         = [];
cfg.elec    = elec;
cfg.channel = 33:42; % 10 EEG channels
layout      = ft_prepare_layout(cfg);

cfg        = [];
cfg.layout = layout;
ft_layoutplot(cfg)

% Preprocess EEG
%--------------------------------------------------------------------------

% Number of bipolar channels
combinations  = nchoosek(33:37,2);
combinations2 = cell(size(combinations));

for c=1:numel(combinations)
    combinations2{c} = ['EEG0',num2str(combinations(c))];
end

C    = size(combinations,1);
data = cell(1,3);
trl  = cell(1,3);
avgs = cell(C,3);

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
    cfg.channel      = 36:45; % 10 EEG channels
    % cfg.channel      = 4:35; % 32 EEG channels
    cfg.continuous   = 'yes';
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = [16,120];
    cfg.dftfilter    = 'yes';       % enable notch filtering to eliminate power line noise
    cfg.dftfreq      = [50,100];    % set up the frequencies for notch filtering
    cfg.dftreplace   = 'zero';      % 'zero' is default'
    cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
    data{cidx}       = ft_preprocessing(cfg); 
end

for c = 1:C % loop over bipolar channels
    for cidx = 1:3 % loop over conditions

        % rereference
        %------------
        cfg = [];
        cfg.channel    = combinations2(c,:);
        cfg.reref      = 'yes';
        cfg.refchannel = combinations2{c,2};
        data1          = ft_preprocessing(cfg,data{cidx});
    
        cfg         = [];
        cfg.channel = combinations2{c,1};
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
        cfg          = [];
        avgs{c,cidx} = ft_timelockanalysis(cfg, epoched1);
    
        clear data1 epoched1
    end
end

% Visualization - single channel
%--------------------------------------------------------------------------
figure
for c=1:C
    cfg           = [];
    cfg.figure    = subplot(2,5,c);
    cfg.linecolor = 'rbk';
    ft_singleplotER(cfg, avgs{c,:});
    legend({conditions{1};conditions{2};conditions{3}});
    title([combinations2{c,1},'-',combinations2{c,2}])
end

mini = 1;
maxi = 0;
for c=1:C
    for cidx=1:3
        a = min(avgs{c,cidx}.avg);
        b = max(avgs{c,cidx}.avg);
        if a<mini
            mini = a;
        end
        if b>maxi
            maxi = b;
        end

    end
end

time2plot = [-10,150];
axisvec = horzcat(time2plot,[mini,maxi]);
color   = {'r','b','k'};
figure
for c=1:C
    subplot(2,5,c);
    for cidx = 1:3
        hold on
        plot(avgs{c,cidx}.time*1000,avgs{c,cidx}.avg,color{cidx});
    end
    xlabel('t/ms')
    
    legend({conditions{1};conditions{2};conditions{3}});
    axis(axisvec) 
    title([combinations2{c,1},'-',combinations2{c,2}])
end


%% Cleaning
rmpath('..')
