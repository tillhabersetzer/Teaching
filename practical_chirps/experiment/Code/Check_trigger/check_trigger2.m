close all
clear all
clc

% path_dataset = fullfile('M:','MEG_Lab',Data','check_trigger','211105','testrun.fif');
% path_dataset = fullfile('M:','MEG_Lab',Data','check_trigger','click.fif');
% path_dataset = fullfile('M:','MEG_Lab',Data','check_trigger','click_loadandwait.fif');
% path_dataset = fullfile('M:','MEG_Lab',Data','check_trigger','upchirps_v2.fif');
path_dataset = fullfile('M:','MEG_Lab','Analysis','Data','bids','sourcedata','subject01','meg','upchirps.fif');

%% check rank of data
% should be less than 306 if maxfiltered

cfg         = [];
cfg.dataset = path_dataset;
cfg.channel = 'meg';
data_select = ft_preprocessing(cfg);
r           = rank(data_select.trial{1}*data_select.trial{1}');


%% check trigger sequence
hdr   = ft_read_header(path_dataset, 'checkmaxfilter', 'no');
event = ft_read_event(path_dataset, 'checkmaxfilter', 'no');

smp = [event.sample];
typ = {event.type};
val = [event.value];

% plot events in scatter plot (all types)
figure
plot(smp/hdr.Fs, val, '.');
title('all event values')

% plot events in scatter plot (only sum channel)
idx            = strcmp(typ,'STI101');
val_sumchannel = val(idx);
figure
subplot(2,1,1)
plot(smp(idx)/hdr.Fs,val_sumchannel,'x')
title(['STI101 events found: ' num2str(sum(idx))])
grid on
subplot(2,1,2)
plot(unique(val(idx)),'x') % already sorted
title(['sorted unique values :' num2str(length(unique(val(idx))))])
grid on

% trigger table
sample  = unique(smp)';
latency = (sample-1)/hdr.Fs;
type    = unique(typ)';

trigarray = nan(length(sample), length(type));

for i=1:numel(sample)
  sel = find(smp==sample(i));
  for j=1:numel(sel)
      trigarray(i, strcmp(type, typ{sel(j)})) = val(sel(j));
  end
end

trigtable = array2table(trigarray, 'VariableNames', type);
trigtable = [table(sample, latency) trigtable];

%% check jitter
jitter = diff(sample)*1000/hdr.Fs; % in ms
histogram(jitter)

%% check trigger alignment
cfg         = [];
cfg.dataset = path_dataset;
data        = ft_preprocessing(cfg);

idx1 = find(contains(data.label,{'MISC'}));
idx2 = find(contains(data.label,{'STI101'}));

figure
subplot(3,1,1)
plot(data.time{1},data.trial{1}(idx1(1),:))
title('MISC Left Audio')
subplot(3,1,2)
plot(data.time{1},data.trial{1}(idx1(2),:))
title('MISC Right Audio')
subplot(3,1,3)
% only plot nonzero elements
idx3 = find(data.trial{1}(idx2,:)~=0);
plot(data.time{1}(idx3),data.trial{1}(idx2,idx3),'x')
title('STI101')

% adapted plot for plotting trigger and misc together

% scale data
a = max(abs(data.trial{1}(idx1(1),:))); % MISC
b = max(abs(data.trial{1}(idx2,idx3))); % STI101
figure('name','scaled data')
plot(data.time{1},1/a*data.trial{1}(idx1(1),:))
hold on
plot(data.time{1}(idx3),1/b*data.trial{1}(idx2,idx3),'rx')
legend({data.label{idx1(1)},data.label{idx2}})

% plot data on top of each other
figure('name','scaled data')
plot(data.time{1},1/a*data.trial{1}(idx1(1),:))
hold on
plot(data.time{1}(idx3),zeros(1,length(data.time{1}(idx3))),'rx')
legend({data.label{idx1(1)},data.label{idx2}})
%xlim([0,10])

