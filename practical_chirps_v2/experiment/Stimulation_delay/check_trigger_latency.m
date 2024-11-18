close all;
clear all;
clc;

% If the file is bigger than 2Gb, it is splitted in two parts

path_dataset1 = fullfile('M:','MEG_Lab','Analysis','Data','stimulation_delay','211125','clicks.fif');
path_dataset2 = fullfile('M:','MEG_Lab','Analysis','Data','stimulation_delay','211125','clicks-1.fif');

%% Concatenate Files
cfg                = [];
cfg.dataset        = path_dataset1;
cfg.checkmaxfilter = false;
cfg.channel        = {'all','-meg'};
data1              = ft_preprocessing(cfg);
cfg.dataset        = path_dataset2;
data2              = ft_preprocessing(cfg);

cfg           =  [];
data          = ft_appenddata(cfg,data1,data2);
trial         = horzcat(data.trial{1},data.trial{2});
data.trial    = [];
data.trial{1} = trial;

% Modify time vector
L             = length(data.time{2});
time2add      = data.time{1}(end)+1/5000:1/5000:data.time{1}(end)+1/5000*L;
time          = horzcat(data.time{1},time2add);
data.time     = [];
data.time{1}  = time;
clear data1 data2 trial time time2add L

%% Number of events

idx_tr = find(contains(data.label,{'STI101'}));  % trigger

events = find(data.trial{1}(idx_tr,1:end-1)-data.trial{1}(idx_tr,2:end)>0);
jitter = diff(events/5); % in ms
histogram(jitter)

%% Plot Channels

idx_al = find(contains(data.label,{'MISC001'})); % original audio left
idx_ra = find(contains(data.label,{'MISC004'})); % recorded audio with latency

% Plot separately
%----------------
figure
subplot(3,1,1)
plot(data.time{1},data.trial{1}(idx_al,:))
title('MISC Left Audio')
xlim([100,105])
subplot(3,1,2)
plot(data.time{1},data.trial{1}(idx_ra,:))
title('MISC Recorded Audio')
xlim([100,105])
subplot(3,1,3)
% only plot nonzero elements
idx_nt = find(data.trial{1}(idx_tr,:)~=0); % no trigger
plot(data.time{1}(idx_nt),data.trial{1}(idx_tr,idx_nt),'x')
title('STI101')
xlim([100,105])

% Plot together
%--------------
% scale data
s_al = max(abs(data.trial{1}(idx_al,:))); 
s_ra = max(abs(data.trial{1}(idx_ra,:))); 
% s_tr = max(abs(data.trial{1}(idx_tr,:))); 

% in ms
figure
hold on
plot(data.time{1}*1000,data.trial{1}(idx_al,:)/s_al,'b')
plot(data.time{1}*1000,data.trial{1}(idx_ra,:)/s_ra,'r')
% only plot nonzero elements
idx_nt = find(data.trial{1}(idx_tr,:)~=0); % no trigger
plot(data.time{1}(idx_nt)*1000,zeros(1,length(idx_nt)),'kx')
title('STI101')
legend({'Audio','Recorded Audio','Trigger'})
% xlim([100,101])
xlim([20*1000,21*1000])


%% Find onsets

% Recorded Audio
%---------------

% normalize (z-value)
% recaudio = data.trial{1}(idx_ra,:);
% recaudio = (recaudio-mean(recaudio))/std(recaudio);

% Or use Initial Audio instead 
%-----------------------------
recaudio = data.trial{1}(idx_al,:);
recaudio = (recaudio-mean(recaudio))/std(recaudio);


% threshold signal 
L       = length(recaudio);
counter = 1;
while counter <= L

    if abs(recaudio(counter))>0.2 % abs(recaudio(counter))>0.5
        recaudio(counter) = 1;
        duration = 0.1*5000; % jump 100 ms, ignore wiggly signal part
        recaudio(counter+1:counter+duration) = 0;
        counter = counter + duration;
    else
        recaudio(counter) = 0;
        counter           = counter +1;
    end

end

onsets_audio = find(recaudio>0);

% figure
% plot(recaudio)

% Trigger
%--------
trigger        = data.trial{1}(idx_tr,:);
onsets_trigger = find(trigger(1:end)-[0,trigger(1:end-1)]>0);

% check onsets: should be the same
disp(['number of events: ',num2str(length(events))])
disp(['number of trigger onsets: ', num2str(length(onsets_trigger))])
disp(['number of audio onsets: ', num2str(length(onsets_audio))])

clear recaudio trigger

%% Latencies

% plot onsets
%------------
% Plot separately
%----------------
figure
subplot(2,1,1)
hold on
plot(data.time{1},data.trial{1}(idx_ra,:),'k')
plot(data.time{1}(onsets_audio),data.trial{1}(idx_ra,onsets_audio),'ro')
title('MISC Recorded Audio')
xlim([100,101])
subplot(2,1,2)
hold on
% only plot nonzero elements
idx_nt = find(data.trial{1}(idx_tr,:)~=0); % no trigger
plot(data.time{1}(idx_nt),data.trial{1}(idx_tr,idx_nt),'kx')
plot(data.time{1}(onsets_trigger),data.trial{1}(idx_tr,onsets_trigger),'ro')
title('STI101')
xlim([100,101])

latency = (onsets_audio - onsets_trigger)*1000/5000; % in ms
figure
histogram(latency)
