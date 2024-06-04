close all
clear 
clc 

% Computations of MMN results for MEG and cEEGrid data
% Analyze impact of trial rejection based on MEG data on MaxFilter

% checkout: https://www.fieldtriptoolbox.org/workshop/natmeg/preprocessing/

% Results
% Application of maxfilter results in significant differences only for
% magnetometers. In that case maxfilter helps to reduce artifacts.
% The gradiometers look almost the same.
% Trial rejection based on MEG channels helps to improve MMN signals

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

% acoustic latency correction for triggers
latency_correction = settings.latency_correction;

%% Script settings
%--------------------------------------------------------------------------
subidx  = 2;
subject = ['sub-',num2str(subidx,'%02d')];

settings.maxfilter = false;
%--------------------------------------------------------------------------

% Path for data
%------------------------------------------------------------------------------------------------------------------------------------
if settings.maxfilter
    datapath = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-mmn-raw_tsss.fif']);
else
    datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-mmn.fif']);
end

%% data overview
%--------------------------------------------------------------------------
% cfg            = [];
% cfg.dataset    = datapath;
% cfg.continuous = 'yes';
% cfg.channel    = 'megmag';
% cfg.viewmode   = 'vertical';
% cfg.blocksize  = 1;                             
% ft_databrowser(cfg);
% 
% cfg            = [];
% cfg.dataset    = datapath;
% cfg.continuous = 'yes';
% cfg.channel    = 'megplanar';
% cfg.viewmode   = 'vertical';
% cfg.blocksize  = 1;                             
% ft_databrowser(cfg);
% 
% cfg                = [];
% cfg.dataset        = datapath;
% cfg.continuous     = 'yes';
% cfg.channel        = 'EEG';
% cfg.viewmode       = 'vertical';
% cfg.blocksize      = 1;   
% cfg.preproc.demean = 'yes';    
% ft_databrowser(cfg);

%% Computation of evoked fields (MEG)
%--------------------------------------------------------------------------

conditions = {'standard','deviant','difference'};
avgs_meg   = cell(1,3);

cfg              = [];
cfg.dataset      = datapath;
cfg.channel      = 'meg'; 
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

% Epoch data
%-----------
cfg     = [];
cfg.trl = trl;  
epoched = ft_redefinetrial(cfg,data); 

trls    = cell(1,2);
trls{1} = find(epoched.trialinfo(:,1) == 1);
trls{2} = find(epoched.trialinfo(:,1) == 2);

% Demean trials
%--------------
cfg                = [];
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.2 0];
epoched            = ft_preprocessing(cfg,epoched);  

% Semi automatic artifact rejection
%----------------------------------
% separately for magnetometers
cfg             = [];
cfg.metric      = 'zvalue';
cfg.layout      = 'neuromag306all.lay';
cfg.channel     = 'megmag';
cfg.keepchannel = 'yes';  % This keeps those channels that are not displayed in the data
epoched         = ft_rejectvisual(cfg,epoched);

% separately for gradiometers
cfg.channel = 'megplanar';
epoched     = ft_rejectvisual(cfg,epoched);

% Downsample data
%----------------
% cfg            = [];
% cfg.resamplefs = 100;
% epoched        = ft_resampledata(cfg,epoched); 

% Timelockanalysis
%-----------------
trls_new    = cell(1,2);
trls_new{1} = find(epoched.trialinfo(:,1) == 1);
trls_new{2} = find(epoched.trialinfo(:,1) == 2);
trialinfo   = [length(trls_new{1}),length(trls_new{2})];

% standard tones
cfg         = [];
cfg.trials  = trls_new{1};
avgs_meg{1} = ft_timelockanalysis(cfg, epoched);

% deviant tones
cfg         = [];
cfg.trials  = trls_new{2};
avgs_meg{2} = ft_timelockanalysis(cfg, epoched);

clear epoched

% differences
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
cfg.operation = 'x1 - x2'; % see https://www.fieldtriptoolbox.org/workshop/natmeg/dipolefitting/
avgs_meg{3}   = ft_math(cfg, avgs_meg{2}, avgs_meg{1});

% Add latency correction 
%-----------------------
% ~5ms 
% 1.8ms transmission delay of sound + max 3ms Recording of analog MEG signals
% -> digital triggers arrive 5 ms to early in MEG recordings
for n = 1:3
    avgs_meg{n}.time  = avgs_meg{n}.time - latency_correction; % in sec
end

%% Visualize data - timeseries
% gradiometers look way superior

% layout = 'neuromag306mag.lay';
layout = 'neuromag306planar.lay';

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbk';
ft_multiplotER(cfg, avgs_meg{1}, avgs_meg{2}, avgs_meg{3});
legend(conditions);

%% Visualize data - topography

layout = 'neuromag306mag.lay';
% layout = 'neuromag306planar.lay';

cfg                 = [];
cfg.layout          = layout;
cfg.xlim            = [0.08 0.15];
cfg.style           = 'straight';
cfg.comment         = 'no';
cfg.marker          = 'off';
cfg.colorbar        = 'southoutside';

cfg.figure = subplot(1,3,1);
ft_topoplotER(cfg, avgs_meg{1});
title('standard');

cfg.figure = subplot(1,3,2);
ft_topoplotER(cfg, avgs_meg{2});
title('deviant');

 cfg.figure= subplot(1,3,3);
ft_topoplotER(cfg, avgs_meg{3});
title('difference');

%% Combining planar gradiometers
% To help with identifying underlying sources we should make use of the 
% other channels in the data. The planar gradiometers are often more easily 
% interpreted, because they are most sensitive right above a source. 

% Combine planar
avgs_cmb = cell(1,3);
cfg     = [];
for n = 1:3
    avgs_cmb{n} = ft_combineplanar(cfg, avgs_meg{n});
end

% Plot results
%-------------
cfg        = [];
cfg.layout = 'neuromag306cmb.lay';
ft_multiplotER(cfg, avgs_cmb{1}, avgs_cmb{2}, avgs_cmb{3});
legend(conditions);

% topographic distribution
cfg                 = [];
cfg.layout          = 'neuromag306cmb.lay';
cfg.xlim            = [0.08 0.15];
cfg.style           = 'straight';
cfg.comment         = 'no';
cfg.marker          = 'off';
cfg.colorbar        = 'southoutside';

cfg.figure = subplot(1,3,1);
ft_topoplotER(cfg, avgs_cmb{1});
title(conditions{1});

cfg.figure = subplot(1,3,2);
ft_topoplotER(cfg, avgs_cmb{2});
title(conditions{2});

cfg.figure= subplot(1,3,3);
ft_topoplotER(cfg, avgs_cmb{3});
title(conditions{3});

%% Computation of evoked potentials (cEEGrid)
%--------------------------------------------------------------------------

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
avgs_left    = cell(2,C(1),3);
avgs_right   = cell(2,C(1),3);
avgs_between = cell(2,C(2),3);

% Computations
%--------------------------------------------------------------------------

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

for aidx = 1:2 % loop over artifact rejections

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
    
            switch aidx
                case 1 % all trials
                    keep_trials = trls;
                case 2 % with trial rejection
                    keep_trials = trls_new;
            end

            % Timelockanalysis
            %-----------------
            % standard tones
            cfg        = [];
            cfg.trials = keep_trials{1};
            avg_sta    = ft_timelockanalysis(cfg, epoched1);
            
            % deviant tones
            cfg        = [];
            cfg.trials = keep_trials{2};
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
                avgs_right{aidx,c,1} = avg_sta;
                avgs_right{aidx,c,2} = avg_dev;
                avgs_right{aidx,c,3} = avg_diff;
            elseif eidx==2
                avgs_left{aidx,c,1} = avg_sta;
                avgs_left{aidx,c,2} = avg_dev;
                avgs_left{aidx,c,3} = avg_diff;
            elseif eidx==3
                avgs_between{aidx,c,1} = avg_sta;
                avgs_between{aidx,c,2} = avg_dev;
                avgs_between{aidx,c,3} = avg_diff;
            end 
    
            clear data1 epoched1 avg_sta avg_dev avg_diff
    
        end
    end
end

%% Visualize cEEGrid results
%--------------------------------------------------------------------------

artifact_correction = {'no','yes'};
% limits
mini = 1;
maxi = 0;
for aidx=1:2
    for c=1:10
        for cidx=1:3
            a = min([avgs_left{aidx,c,cidx}.avg,avgs_right{aidx,c,cidx}.avg]);
            b = max([avgs_left{aidx,c,cidx}.avg,avgs_right{aidx,c,cidx}.avg]);
            if a<mini
                mini = a;
            end
            if b>maxi
                maxi = b;
            end
        end
    end
end

% add channel combinations between ears
for aidx=1:2
    for c=1:25
        for cidx=1:3
            a = min([avgs_between{aidx,c,cidx}.avg]);
            b = max([avgs_between{aidx,c,cidx}.avg]);
            if a<mini
                mini = a;
            end
            if b>maxi
                maxi = b;
            end
        end
    end
end

time2plot = [-200,800];
axisvec = horzcat(time2plot,[mini,maxi]);
color   = {'r','b','k'};

% right ear
%----------
for aidx=1:2
    figure
    for c=1:10 % loop over channel combinations
        subplot(2,5,c);
        for cidx = 1:3 % loop over conditions
            hold on
            plot(avgs_right{aidx,c,cidx}.time*1000,avgs_right{aidx,c,cidx}.avg,color{cidx});
        end
        xlabel('t/ms')
        
        if c==1
            legend({conditions{1};conditions{2};conditions{3}});
        end
        axis(axisvec) 
        title([combinations_label{1}{c,1},'-',combinations_label{1}{c,2}])
    end
    sgtitle([subject,': right ear (artifact correction: ',artifact_correction{aidx},')'])
end

% left ear
%----------
for aidx=1:2
    figure
    for c=1:10 % loop over channel combinations
        subplot(2,5,c);
        for cidx = 1:3 % loop over conditions
            hold on
            plot(avgs_left{aidx,c,cidx}.time*1000,avgs_left{aidx,c,cidx}.avg,color{cidx});
        end
        xlabel('t/ms')
        
        if c==1
            legend({conditions{1};conditions{2};conditions{3}});
        end
        axis(axisvec) 
        title([combinations_label{2}{c,1},'-',combinations_label{2}{c,2}])
    end
    sgtitle([subject,': left ear (artifact correction: ',artifact_correction{aidx},')'])
end

% between ears
%-------------
for aidx=1:2
    figure
    for c=1:25 % loop over channel combinations
        subplot(5,5,c);
        for cidx = 1:3 % loop over conditions
            hold on
            plot(avgs_between{aidx,c,cidx}.time*1000,avgs_between{aidx,c,cidx}.avg,color{cidx});
        end
        xlabel('t/ms')
        
        if c==1
            legend({conditions{1};conditions{2};conditions{3}});
        end
        axis(axisvec) 
        title([combinations_label{3}{c,1},'-',combinations_label{3}{c,2}])
    end
    sgtitle([subject,': between ears (artifact correction: ',artifact_correction{aidx},')'])
end

%% only visualize the difference between standard and deviant tones
%--------------------------------------------------------------------------

% limits
mini = 1;
maxi = 0;
for aidx=1:2
    for c=1:10       
        a = min([avgs_left{aidx,c,3}.avg,avgs_right{aidx,c,3}.avg]);
        b = max([avgs_left{aidx,c,3}.avg,avgs_right{aidx,c,3}.avg]);
        if a<mini
            mini = a;
        end
        if b>maxi
            maxi = b;
        end      
    end
end

% add channel combinations between ears
for aidx=1:2
    for c=1:25
        a = min([avgs_between{aidx,c,3}.avg]);
        b = max([avgs_between{aidx,c,3}.avg]);
        if a<mini
            mini = a;
        end
        if b>maxi
            maxi = b;
        end
    end
end

time2plot = [-200,800];
axisvec = horzcat(time2plot,[mini,maxi]);
color   = {'r','b'};

% right ear
%----------
figure
for c=1:10 % loop over channel combinations
    for aidx=1:2
        subplot(2,5,c);
        hold on
        plot(avgs_right{aidx,c,3}.time*1000,avgs_right{aidx,c,3}.avg,color{aidx});
        xlabel('t/ms')
        
        if c==1
            legend(artifact_correction);
        end
        axis(axisvec) 
        title([combinations_label{1}{c,1},'-',combinations_label{1}{c,2}])
    end
end
sgtitle([subject,': right ear'])

% left ear
%----------
% right ear
%----------
figure
for c=1:10 % loop over channel combinations
    for aidx=1:2
        subplot(2,5,c);
        hold on
        plot(avgs_left{aidx,c,3}.time*1000,avgs_left{aidx,c,3}.avg,color{aidx});
        xlabel('t/ms')
        
        if c==1
            legend(artifact_correction);
        end
        axis(axisvec) 
        title([combinations_label{2}{c,1},'-',combinations_label{2}{c,2}])
    end
end
sgtitle([subject,': left ear'])


%% Clean-Up
rmpath('..')








