close all
clear 
clc 

% Setup
% sub-01: 42 channel EEG, 32 channel EEG cap + 10 extra electrodes around ear
% sub-02: 10 channel ear EEG

%% Script settings
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

subidx  = 2; 
subject = ['sub-',num2str(subidx,'%02d')];

conditions = {'clicks','upchirps','downchirps'};
%--------------------------------------------------------------------------

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

%% Import data
%--------------------------------------------------------------------------
dir2save   = fullfile(settings.path2project,'derivatives',subject,'sensorlevel');
data       = load(fullfile(dir2save,[subject,'_chirps_erps_earEEG.mat']));
epochs     = data.epochs;
conditions = data.conditions;

clear data

%% Layout ear EEG
%--------------------------------------------------------------------------

switch subject
    case 'sub-01'
        ear_channels = {'EEG038','EEG039','EEG040','EEG041','EEG042','EEG008',... % left
                        'EEG033','EEG034','EEG035','EEG036','EEG037','EEG012'}; % right
        ear_channels_rename = {'L01','L02','L03','L04','L05','L06',...
                               'R01','R02','R03','R04','R05','R06'};
    case 'sub-02'
        ear_channels = {'EEG006','EEG007','EEG008','EEG009','EEG010',... % left
                        'EEG001','EEG002','EEG003','EEG004','EEG005'}; % right
        ear_channels_rename = {'L01','L02','L03','L04','L05',...
                               'R01','R02','R03','R04','R05'};
end

% Plot geomtery
%--------------
datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-clicks.fif']);
elec     = ft_read_sens(datapath, 'senstype', 'eeg');

elec.label = rename_channels(elec.label, ear_channels, ear_channels_rename);
headshape = ft_read_headshape(datapath);

figure
hold on
ft_plot_sens(elec,'label','on','elecshape','disc');
ft_plot_headshape(headshape);

% layout
%-------
cfg         = [];
cfg.elec    = elec;
cfg.channel = 'all'; 
layout      = ft_prepare_layout(cfg);

cfg        = [];
cfg.layout = layout;
ft_layoutplot(cfg)

% Rename EEG channels
%--------------------
for cidx = 1:3
    epochs{cidx}.label = rename_channels(epochs{cidx}.label, ear_channels, ear_channels_rename);
end

%% Computation of bipolar channel configurations / Rereferencing
%--------------------------------------------------------------------------

% Number of channels
n_chan = length(epochs{1}.label);

% Number of possible bipolar channels N*(N-1)/2
n_bip = n_chan*(n_chan-1)/2;

% Compute channel combinations in terms of index
% If channels are missing, the index does not correspond to actual channels
% 3 =~ EEG003 if EEG002 is missing
combinations = zeros(n_chan,2);
counter      = 1;
for cidx1 = 1:n_chan % 1st channel    
    for cidx2 = cidx1+1:n_chan % 2nd channel
        combinations(counter,1) = cidx1;
        combinations(counter,2) = cidx2;
        counter = counter +1;
    end
end

labels             = epochs{2}.label; % Get channel names
combinations_names = cell(n_bip,3);
avgs               = cell(1,3); % 3 conditions 
stds_mean          = cell(1,3);
for conidx = 1:3
    avgs{conidx}      = zeros(n_bip,length(epochs{conidx}.time{1}));
    stds_mean{conidx} = zeros(size(avgs{conidx}));
end

% Get idx for baseline time window
timewin_baseline_idx = find(epochs{1}.time{1} < 0);

% Computation of bipolar ERPs and standard deviations
%----------------------------------------------------
for bidx = 1:n_bip % Loop over bipolar channels

    channels                     = {labels{combinations(bidx,1)},labels{combinations(bidx,2)}};
    combinations_names(bidx,1:2) = channels;
    combinations_names{bidx,3}   = [channels{1},'-',channels{2}];

    for conidx = 1:3 % Loop over conditions
    
        % Subselection of channels
        cfg              = [];
        cfg.channel      = channels;
        epochs_selection = ft_selectdata(cfg,epochs{conidx});
    
        % Rereferencing
        trials = zeros(length(epochs_selection.trial),length(epochs_selection.time{1}));
        for tidx = 1:size(trials,1) % loop over all trials
            trials(tidx,:) = epochs_selection.trial{tidx}(1,:) - epochs_selection.trial{tidx}(2,:);
        end
    
        % Compute average over trials
        avg = mean(trials,1);
        % Computate standard deviation of the mean
        stds_mean{conidx}(bidx,:) = std(trials,[],1)/sqrt(size(trials,1));
    
        % Subtract baseline from average
        avg = avg - mean(avg(timewin_baseline_idx(1):timewin_baseline_idx(end)),2);
    
        avgs{conidx}(bidx,:) = avg;
        clear avg
    
        disp(['Bipolar channel ',num2str(bidx),'/',num2str(n_bip),' computed.'])
    end
end

% Check for polarity and flip avg - Use deviant condition as reference
%---------------------------------------------------------------------
ref_timewin     = [0.01,0.05];
ref_timewin_idx = find(epochs{1}.time{1} >= ref_timewin(1) & epochs{1}.time{1} <= ref_timewin(2));

% Sort based on upchirp amplitde
conidx_reference = find(contains(conditions,'upchirp'));

for bidx = 1:n_bip

    % mean value in certain time perion (should be > 0)
    signal = avgs{conidx_reference}(bidx,ref_timewin_idx(1):ref_timewin_idx(end));
    [~,idx] = max(abs(signal));
    if signal(idx) < 0 % flip 
        avgs{1}(bidx,:) = -1*avgs{1}(bidx,:);
        avgs{2}(bidx,:) = -1*avgs{2}(bidx,:);
        avgs{3}(bidx,:) = -1*avgs{3}(bidx,:);

        % Variance does not change when flipped
        % Var(aX) = a^2*Var(X), a=-1
        % change channel combinations
        channels                     = combinations_names(bidx,1:2);
        combinations_names(bidx,1:2) = {channels{2},channels{1}};
        combinations_names{bidx,3}   = [channels{2},'-',channels{1}];
        combinations(bidx,:)         = fliplr(combinations(bidx,:));
        disp('Reference flipped.')
    end
    clear signal
end

%% Split data into left and right channels
%--------------------------------------------------------------------------
% Left channels (Both channels must be on the left side)
lidx = find(contains(combinations_names(:,1),'L') & contains(combinations_names(:,2),'L'));

% Right channels (Both channels must be on the right side)
ridx = find(contains(combinations_names(:,1),'R') & contains(combinations_names(:,2),'R'));

% Mixed ears
midx = find(contains(combinations_names(:,1),'L') & contains(combinations_names(:,2),'R') | ...
            contains(combinations_names(:,1),'R') & contains(combinations_names(:,2),'L'));

% Check combinations
if ~isequal(length(lidx)+length(ridx)+length(midx),n_bip)
    error('Unexpected number of channel combinations!')
end

% Sort channels
%--------------

% Sort based on upchirp amplitde
cidx = find(contains(conditions,'upchirp'));
vals_left   = max(abs(avgs{cidx}(lidx,ref_timewin_idx(1):ref_timewin_idx(end))),[],2);
vals_right  = max(abs(avgs{cidx}(ridx,ref_timewin_idx(1):ref_timewin_idx(end))),[],2);
vals_linked = max(abs(avgs{cidx}(midx,ref_timewin_idx(1):ref_timewin_idx(end))),[],2);
vals_all    = max(abs(avgs{cidx}(:,ref_timewin_idx(1):ref_timewin_idx(end))),[],2);

[vals_left,idx] = sort(abs(vals_left));
lidx            = lidx(idx);

[vals_right,idx] = sort(abs(vals_right));
ridx             = ridx(idx);

[vals_linked,idx] = sort(abs(vals_linked));
midx              = midx(idx);

[vals_all,idx] = sort(abs(vals_all));
aidx           = 1:n_bip;
aidx           = aidx(idx)';

mini = min(vals_all,[],"all");
maxi = max(vals_all,[],"all");

% Show channels with maximum effect
%----------------------------------
figure
subplot(2,3,1)
plot(1:length(lidx),vals_left,'ro')
xticks(1:length(lidx))
xticklabels(combinations_names(lidx,3));
% xtickangle(45)
% ylabel('U / uV')
ylim([mini,maxi])
title('Left ear')
subplot(2,3,2)
plot(1:length(ridx),vals_right,'go')
xticks(1:length(lidx))
xticklabels(combinations_names(ridx,3));
% xtickangle(45)
ylim([mini,maxi])
title('Right ear')
subplot(2,3,3)
plot(1:length(midx),vals_linked,'bo')
xticks(1:length(midx))
xticklabels(combinations_names(midx,3));
% xtickangle(45)
ylim([mini,maxi])
title('Linked ears')
subplot(2,3,4:6)
% plot(1:length(aidx),vals_all*10^6,'o')
hold on
plot(find(ismember(aidx,lidx)),vals_left,'ro')
plot(find(ismember(aidx,ridx)),vals_right,'go')
plot(find(ismember(aidx,midx)),vals_linked,'bo')
legend({'Left ear','Right ear', 'Linked ears'},'location','northwest')
xticks(1:length(aidx))
xticklabels(combinations_names(aidx,3));
% xtickangle(45)
ylim([mini,maxi])
% Get handle to current axes.
ax = gca;
ax.XAxis.FontSize = 5;

sgtitle("Maximum absolute value")

%% Plot results
%-------------------------------------------------------------------------- 
time_vec = epochs{1}.time{1};
colors   = {'r','b','g'};

x_text = max(time_vec)-0.1;
y_text = {3,2,1};
figure
for conidx = 1:3 % Plot both conditions
    % Left ear
    subplot(3,1,1)
    hold on
    plot(time_vec,avgs{conidx}(lidx,:)*10^6,colors{conidx})
   
    ylabel('U / uV')
    title('Left ear')
    text(x_text, y_text{conidx}, conditions{conidx}, 'Color', colors{conidx})

    % Right ear
    subplot(3,1,2)
    hold on
    plot(time_vec,avgs{conidx}(ridx,:)*10^6,colors{conidx})
    ylabel('U / uV')
    title('Right ear')

    % Linked ears
    subplot(3,1,3)
    hold on
    plot(time_vec,avgs{conidx}(midx,:)*10^6,colors{conidx})
    xlabel('t / s')
    ylabel('U / uV')
    title('Linked ears')
end

%% Other plot 3D - too messy
%--------------------------------------------------------------------------
% plot both conditions with standard errors of mean in 3D

% Choose which channels to plot
ear = 'left';
% ear = 'right';
% ear = 'linked';

switch ear
    case 'left'
        chan_idxs  = lidx;
        title_name = 'Left ear';
    case 'right'
        chan_idxs  = ridx;
        title_name = 'Right ear';
    case 'linked'
        chan_idxs  = midx;
        title_name = 'Linked ears';
end

figure
for conidx = 1:3
    for pidx = 1:length(chan_idxs) % loop over channels
    
        chan_idx = chan_idxs(pidx); % channel index
        hold on
        plot3(pidx*ones(size(time_vec)),time_vec,avgs{conidx}(chan_idx,:)*10^6,colors{conidx},'linewidth',1)
        % add standard deviation of mean
        patch([pidx*ones(size(time_vec)), pidx*ones(size(time_vec))],[time_vec fliplr(time_vec)], [(avgs{conidx}(chan_idx,:)-stds_mean{conidx}(chan_idx,:))*10^6 fliplr((avgs{conidx}(chan_idx,:)+stds_mean{conidx}(chan_idx,:)))*10^6], colors{conidx}, 'FaceAlpha',0.1, 'EdgeColor','none','HandleVisibility','off')
    
    end
end
xticks(1:length(chan_idxs))
xticklabels(combinations_names(chan_idxs,3));
xlabel('bipolar channels')
ylabel('t / s')
zlabel('U / uV')
grid on
title(title_name)
view(3)

%% Plot channels with largest effect 
%--------------------------------------------------------------------------

figure
for i = 1:12
    subplot(3,4,i) % both conditions
    chan_idx = aidx(end+1-i);
    label    = combinations_names(chan_idx,3);
    for conidx = 1:3
        hold on
        plot(time_vec,avgs{conidx}(chan_idx,:)*10^6,colors{conidx})
        patch([time_vec fliplr(time_vec)], [(avgs{conidx}(chan_idx,:)-stds_mean{conidx}(chan_idx,:))*10^6 fliplr((avgs{conidx}(chan_idx,:)+stds_mean{conidx}(chan_idx,:)))*10^6], colors{conidx}, 'FaceAlpha',0.1, 'EdgeColor','none','HandleVisibility','off')
    end
    xlabel('t / s')
    ylabel('U / uV')
    legend(conditions)
    title([label," (Max abs:",num2str(vals_all(end+1-i)),')'])
end  

%% Plot single channel - Select channels based on effect size
%--------------------------------------------------------------------------
chan_idx = aidx(end);
label    = combinations_names(chan_idx,3);
% label    = {'R06-L03'};
chan_idx = find(contains(combinations_names(:,3),label));

figure
for conidx = 1:3
    hold on
    plot(time_vec,avgs{conidx}(chan_idx,:)*10^6,colors{conidx})
    patch([time_vec fliplr(time_vec)], [(avgs{conidx}(chan_idx,:)-stds_mean{conidx}(chan_idx,:))*10^6 fliplr((avgs{conidx}(chan_idx,:)+stds_mean{conidx}(chan_idx,:)))*10^6], colors{conidx}, 'FaceAlpha',0.1, 'EdgeColor','none','HandleVisibility','off')
end
xlabel('t / s')
ylabel('U / uV')
legend(conditions)
title(['Both conditions, single channel (',label,')'])

%% Cleaning
rmpath('..')

%% functions

function data_out = rename_channels(data_in,old_labels,new_labels)
  
% data_out: ioutput label vector with adjusted channel names
% data_in: input label vector 
% old_labels, new_labels: label names that should be adjusted

    data_out = data_in;

    % Rename channels
    %----------------   
    C = length(old_labels);
    idx = zeros(1,C);
    for i = 1:C
        idx(i) = find(contains(data_in,old_labels{i}));
    end
    for i=1:C
        data_out{idx(i)} = new_labels{i};
    end

end
