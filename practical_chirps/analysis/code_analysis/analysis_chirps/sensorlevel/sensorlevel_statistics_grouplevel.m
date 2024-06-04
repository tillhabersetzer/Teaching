close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
for i = 3:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjects{i-2} = [subject,num2str(i)]; 
end

% remove bad subjects
subjects(8) = []; subjects(17) = []; subjects(18) = [];

check = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'));

% load preprocessed data
%-----------------------
N_subj       = length(subjects);
avg_pre_mag  = cell(1,N_subj);
avg_pre_grad = cell(1,N_subj);
avg_pre_cmb  = cell(1,N_subj);
avg_pst_mag  = cell(1,N_subj);
avg_pst_grad = cell(1,N_subj);
avg_pst_cmb  = cell(1,N_subj);

info    = cell(1,N_subj);
trials  = zeros(1,N_subj);

for s = 1:N_subj
    eval(subjects{s})
    data            = importdata(fullfile(subjectdata.chirps_sensorlevel,[subjectdata.subjectname,'_averages.mat']));
    avg             = data.avg;
    % normalize data
    %---------------
    selmag  = ft_chantype(avg.label, 'megmag');
    selgrad = ft_chantype(avg.label, 'megplanar');
    avg.avg(selmag,:) = avg.avg(selmag,:)./max(avg.avg(selmag,:),[],'all');
    avg.avg(selgrad,:) = avg.avg(selgrad,:)./max(avg.avg(selgrad,:),[],'all');
    
    cfg             = [];
    cfg.latency     = [-0.45 0];
    cfg.channel     = 'megmag';
    avg_pre_mag{s}  = ft_selectdata(cfg,avg);
    cfg.channel     = 'megplanar';
    avg_pre_grad{s} = ft_selectdata(cfg,avg);
    avg_pre_cmb{s}  = ft_selectdata(cfg,avg);
    
    cfg             = [];
    cfg.latency     = [0 0.45];
    cfg.channel     = 'megmag';
    avg_pst_mag{s}  = ft_selectdata(cfg,avg);
    cfg.channel     = 'megplanar';
    avg_pst_grad{s} = ft_selectdata(cfg,avg);
    avg_pst_cmb{s}  = ft_selectdata(cfg,avg);
    
    cfg                = [];
    cfg.channel        = 'megplanar';
    cfg.method         = 'sum';
    cfg.updatesens     = 'yes'; 
    cfg.demean         = 'yes';
    cfg.baselinewindow = [-0.5,0]; % only demean by baselinewindow to do offset correction 
    avg_cmb            = ft_combineplanar(cfg,avg);
    cfg                = [];
    cfg.latency        = [-0.45 0];
    avg_pre_cmb{s}     = ft_selectdata(cfg,avg_cmb);
    cfg.latency        = [0 0.45];
    avg_pst_cmb{s}     = ft_selectdata(cfg,avg_cmb);
    
    avg_pre_mag{s}.time  = avg_pst_mag{s}.time; % change time axis
    avg_pre_grad{s}.time = avg_pst_mag{s}.time; % change time axis
    avg_pre_cmb{s}.time  = avg_pst_mag{s}.time; % change time axis

    clear avg avg_cmb
    info{s}   = data.info;
    trials(s) = info{s}.number_of_epochs;
    
    % optional - use abs values
    %--------------------------
%     cfg           = [];
%     cfg.operation = 'abs';
%     cfg.parameter = 'avg';
%     avg_pre{s}    = ft_math(cfg,avg_pre{s});
%     avg_pst{s}    = ft_math(cfg,avg_pst{s});
    
    % optional - test against averaged baseline -> looks better without
    %------------------------------------------
    % average over time (columns) and replace matrix data with average
%      avg_pre{s}.avg = repmat(mean(avg_pre{s}.avg,2),1,size(avg_pre{s}.avg,2));
end

% calculate grandaverage 
%-----------------------
cfg           = [];
cfg.latency   = 'all';
gavg_pre_mag  = ft_timelockgrandaverage(cfg,avg_pre_mag{:});
gavg_pre_grad = ft_timelockgrandaverage(cfg,avg_pre_grad{:});
gavg_pre_cmb  = ft_timelockgrandaverage(cfg,avg_pre_cmb{:});
gavg_pst_mag  = ft_timelockgrandaverage(cfg,avg_pst_mag{:});
gavg_pst_grad = ft_timelockgrandaverage(cfg,avg_pst_grad{:});
gavg_pst_cmb  = ft_timelockgrandaverage(cfg,avg_pst_cmb{:});

cfg            = [];
cfg.operation  = 'subtract';
cfg.parameter  = 'avg';
raweffect_mag  = ft_math(cfg,gavg_pst_mag,gavg_pre_mag);
raweffect_grad = ft_math(cfg,gavg_pst_grad,gavg_pre_grad);
raweffect_cmb  = ft_math(cfg,gavg_pst_cmb,gavg_pre_cmb);

if check
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306mag_helmet.mat';
    ft_multiplotER(cfg,gavg_pre_mag,gavg_pst_mag,raweffect_mag);
    sgtitle('grandaverage: magnetometers') 
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306planar_helmet.mat';
    ft_multiplotER(cfg,gavg_pre_grad,gavg_pst_grad,raweffect_grad);
    sgtitle('grandaverage: gradiometers') 
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb_helmet.mat';
    ft_multiplotER(cfg,gavg_pre_cmb,gavg_pst_cmb,raweffect_cmb);
    sgtitle('grandaverage: combined gradiometers') 
    
    % to get an overview about subject level data
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb_helmet.mat';
    ft_multiplotER(cfg,avg_pst_cmb{:});
    sgtitle('all subjects: combined gradiometers') 
    
    n = 6;
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb_helmet.mat';
    ft_multiplotER(cfg,avg_pre_cmb{n},avg_pst_cmb{n});
    sgtitle([subjects{n},': combined gradiometers']) 
end

%% statistic
%-----------
cfg                         = [];
cfg.method                  = 'montecarlo';
cfg.statistic               = 'ft_statfun_depsamplesT';
cfg.correctm                = 'cluster';
cfg.clusteralpha            = 0.05;
cfg.clusterstatistic        = 'maxsum';
cfg.minnbchan               = 3;
cfg.tail                    = 0;
cfg.clustertail             = 0;
cfg.alpha                   = 0.025;
cfg.numrandomization        = 1000;
design                      = zeros(2,2*N_subj);
design(1,1:N_subj)          = 1;
design(1,N_subj+1:2*N_subj) = 2;
design(2,1:N_subj)          = [1:N_subj];
design(2,N_subj+1:2*N_subj) = [1:N_subj];
cfg.design                  = design;
cfg.ivar                    = 1; % row of design matrix that contains independent variable 
cfg.uvar                    = 2; % row of design matrix that contains unit of observation

cfg.channel         = 'megmag';
cfg_neighb          = [];
cfg_neighb.channel  = 'megmag';
cfg_neighb.method   = 'template';
cfg_neighb.template = 'neuromag306mag_neighb.mat';
cfg.neighbours      = ft_prepare_neighbours(cfg_neighb,avg_pre_mag{1}.grad);
%ft_neighbourplot(cfg,avg_pre_mag{1})
stat_mag            = ft_timelockstatistics(cfg,avg_pst_mag{:},avg_pre_mag{:});

cfg.channel         = 'megplanar';
cfg_neighb          = [];
cfg_neighb.channel  = 'megplanar';
cfg_neighb.method   = 'template';
cfg_neighb.template = 'neuromag306planar_neighb.mat';
cfg.neighbours      = ft_prepare_neighbours(cfg_neighb,avg_pre_grad{1}.grad);
%ft_neighbourplot(cfg,avg_pre_grad{1})
stat_grad           = ft_timelockstatistics(cfg,avg_pst_grad{:},avg_pre_grad{:});

% change to one sided test for combined gradiometers - only positive
% clusters
cfg.tail        = 1;
cfg.clustertail = 1;
cfg.alpha       = 0.05;

cfg.channel         = 'megplanar';
cfg_neighb          = [];
cfg_neighb.channel  = 'megplanar';
cfg_neighb.method   = 'template';
cfg_neighb.template = 'neuromag306cmb_neighb.mat';
cfg.neighbours      = ft_prepare_neighbours(cfg_neighb,avg_pre_cmb{1}.grad);
%ft_neighbourplot(cfg,avg_pre_cmb{1})
stat_cmb            = ft_timelockstatistics(cfg,avg_pst_cmb{:},avg_pre_cmb{:});

% directory to save graphics
%---------------------------
new_dir = fullfile('Z:','analysis','generated_data','grouplevel','chirps','sensorlevel');
if ~exist(new_dir, 'dir')
   mkdir(new_dir)
end

% time windows for topoplots
windows = [0,0.05;0.05,0.075;0.075,0.1;0.1,0.125;0.125,0.15;0.15,0.175;0.175,0.2;...
           0.2,0.225;0.225,0.25;0.25,0.3;0.3,0.35;0.35,0.4;0.4,0.45];

%% make nice plots - magnetometers
%---------------------------------       
pos_cluster_pvals = [stat_mag.posclusters(:).prob];
pos_signif_clust  = find(pos_cluster_pvals < stat_mag.cfg.alpha); 
pos               = ismember(stat_mag.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals = [stat_mag.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < stat_mag.cfg.alpha);
neg               = ismember(stat_mag.negclusterslabelmat, neg_signif_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect_mag.label, stat_mag.label);
% find correct indices for timewindows
idx     = dsearchn(stat_mag.time',reshape(windows,1,numel(windows))');
idx     = reshape(idx,size(windows,1),size(windows,2));
    
% Plot over raweffect
%--------------------
figure('Position', get(0, 'Screensize'))
for k = 1:13
   subplot(3,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   pos_int  = zeros(numel(raweffect_mag.label),1);
   neg_int  = zeros(numel(raweffect_mag.label),1);
   %pos_int(i1) = all(pos(i2,idx(k,1):idx(k,2)), 2);
   %neg_int(i1) = all(neg(i2,idx(k,1):idx(k,2)), 2);
   % less severe
   pos_int(i1) = any(pos(i2,idx(k,1):idx(k,2)), 2);
   neg_int(i1) = any(neg(i2,idx(k,1):idx(k,2)), 2);
   cfg.highlight        = 'on';
   % Get the index of each significant channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 10;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306mag_helmet.mat';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg,raweffect_mag);
   set(gca,'fontsize', 12)
end
sgtitle('grandaverage - raweffect magnetometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_topoplot_timeseries_raweffect_magnetometers.png'))

% Plot over t-values
%-------------------
% change t-values to abs for plotting
% stati_mag      = stat_mag;
% stati_mag.stat = abs(stat_mag.stat);

figure('Position', get(0, 'Screensize'))
for k = 1:13
   subplot(3,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   cfg.zlim = [-5,5];
   %cfg.zlim = [0,5];
   %cfg.zlim = [0 max(stati.stat,[],'all')];
   cfg.highlight        = 'on';
   % int = all(stat_mag.mask(:,idx(k,1):idx(k,2)), 2);
   % less severe
   int = any(stat_mag.mask(:,idx(k,1):idx(k,2)), 2); 
   cfg.highlightchannel = find(int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 12;
   cfg.colormap         = parula;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306mag_helmet.mat';
   cfg.interactive      = 'no';
   cfg.parameter        = 'stat';
   ft_topoplotER(cfg, stat_mag);
   set(gca,'ColorScale','linear') % if you want to plot on a logarithmic scale
   set(gca,'fontsize', 12)
end
posi        = get(subplot(3,5,13),'Position'); %[left bottom width height]
posi(1)     = posi(1)+posi(3)+0.01; posi(3) = 0.25*posi(3);
c           = colorbar('Position',posi);
c.LineWidth = 1;
c.FontSize  = 15;
title(c,'t-val','fontweight','bold');
%set(c,'Ticks',[0.1,1,2,3,4,5],'TickLabels',{'0.1','1','2','3','4','5'})

sgtitle('grandaverage - t-values - magnetometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_topoplot_timeseries_tvalues_magnetometers.png'))

% timeseries
cfg               = [];
cfg.layout        = 'neuromag306mag_helmet.mat';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor    = 'r';
figure('Position', get(0, 'Screensize'))
ft_multiplotER(cfg,stat_mag);

sgtitle('grandaverage - timeseries magnetometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_timeseries_magnetometers.png'))


%% make nice plots - gradiometers
%---------------------------------       
pos_cluster_pvals = [stat_grad.posclusters(:).prob];
pos_signif_clust  = find(pos_cluster_pvals < stat_grad.cfg.alpha); 
pos               = ismember(stat_grad.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals = [stat_grad.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < stat_grad.cfg.alpha);
neg               = ismember(stat_grad.negclusterslabelmat, neg_signif_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect_grad.label, stat_grad.label);
% find correct indices for timewindows
idx     = dsearchn(stat_grad.time',reshape(windows,1,numel(windows))');
idx     = reshape(idx,size(windows,1),size(windows,2));
    
% Plot over raweffect
%--------------------
figure('Position', get(0, 'Screensize'))
for k = 1:13
   subplot(3,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   pos_int  = zeros(numel(raweffect_grad.label),1);
   neg_int  = zeros(numel(raweffect_grad.label),1);
   %pos_int(i1) = all(pos(i2,idx(k,1):idx(k,2)), 2);
   %neg_int(i1) = all(neg(i2,idx(k,1):idx(k,2)), 2);
   % less severe
   pos_int(i1) = any(pos(i2,idx(k,1):idx(k,2)), 2);
   neg_int(i1) = any(neg(i2,idx(k,1):idx(k,2)), 2);
   cfg.highlight        = 'on';
   % Get the index of each significant channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 10;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306planar_helmet.mat';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg,raweffect_grad);
   set(gca,'fontsize', 12)
end
sgtitle('grandaverage - raweffect gradiometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_topoplot_timeseries_raweffect_gradiometers.png'))

% Plot over t-values
%-------------------
% change t-values to abs for plotting
% stati_grad      = stat_grad;
% stati_grad.stat = abs(stat_grad.stat);

figure('Position', get(0, 'Screensize'))
for k = 1:13
   subplot(3,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   cfg.zlim = [-5,5];
   %cfg.zlim = [0,5];
   %cfg.zlim = [0 max(stati.stat,[],'all')];
   cfg.highlight        = 'on';
   % int = all(stat_grad.mask(:,idx(k,1):idx(k,2)), 2);
   % less severe
   int = any(stat_grad.mask(:,idx(k,1):idx(k,2)), 2); 
   cfg.highlightchannel = find(int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 12;
   cfg.colormap         = parula;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306planar_helmet.mat';
   cfg.interactive      = 'no';
   cfg.parameter        = 'stat';
   ft_topoplotER(cfg, stat_grad);
   set(gca,'ColorScale','linear') % if you want to plot on a logarithmic scale
   set(gca,'fontsize', 12)
end
posi        = get(subplot(3,5,13),'Position'); %[left bottom width height]
posi(1)     = posi(1)+posi(3)+0.01; posi(3) = 0.25*posi(3);
c           = colorbar('Position',posi);
c.LineWidth = 1;
c.FontSize  = 15;
title(c,'t-val','fontweight','bold');
%set(c,'Ticks',[0.1,1,2,3,4,5],'TickLabels',{'0.1','1','2','3','4','5'})

sgtitle('grandaverage - t-values - gradiometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_topoplot_timeseries_tvalues_gradiometers.png'))

% timeseries
cfg               = [];
cfg.layout        = 'neuromag306planar_helmet.mat';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor    = 'r';
figure('Position', get(0, 'Screensize'))
ft_multiplotER(cfg,stat_grad);

sgtitle('grandaverage - timeseries gradiometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_timeseries_gradiometers.png'))


%% make nice plots - combined gradiometers
%-----------------------------------------       
pos_cluster_pvals = [stat_cmb.posclusters(:).prob];
pos_signif_clust  = find(pos_cluster_pvals < stat_cmb.cfg.alpha); 
pos               = ismember(stat_cmb.posclusterslabelmat, pos_signif_clust);

% neg_cluster_pvals = [stat_cmb.negclusters(:).prob];
% neg_signif_clust  = find(neg_cluster_pvals < stat_cmb.cfg.alpha);
% neg               = ismember(stat_cmb.negclusterslabelmat, neg_signif_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect_cmb.label, stat_cmb.label);
% find correct indices for timewindows
idx     = dsearchn(stat_cmb.time',reshape(windows,1,numel(windows))');
idx     = reshape(idx,size(windows,1),size(windows,2));
    
% Plot over raweffect
%--------------------
figure('Position', get(0, 'Screensize'))
for k = 1:13
   subplot(3,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   pos_int  = zeros(numel(raweffect_cmb.label),1);
%    neg_int  = zeros(numel(raweffect_cmb.label),1);
   pos_int(i1) = all(pos(i2,idx(k,1):idx(k,2)), 2);
   %neg_int(i1) = all(neg(i2,idx(k,1):idx(k,2)), 2);
   % less severe
%      pos_int(i1) = any(pos(i2,idx(k,1):idx(k,2)), 2);
%    neg_int(i1) = any(neg(i2,idx(k,1):idx(k,2)), 2);
   cfg.highlight        = 'on';
   % Get the index of each significant channel
%    cfg.highlightchannel = find(pos_int | neg_int);
   cfg.highlightchannel = find(pos_int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 10;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306cmb_helmet.mat';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg,raweffect_cmb);
   set(gca,'fontsize', 12)
end
sgtitle('grandaverage - raweffect combined gradiometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_topoplot_timeseries_raweffect_combined_gradiometers.png'))

% Plot over t-values
%-------------------
% change t-values to abs for plotting
% stati_cmb      = stat_cmb;
% stati_cmb.stat = abs(stat_cmb.stat);

figure('Position', get(0, 'Screensize'))
for k = 1:13
   subplot(3,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   cfg.zlim = [-5,5];
   %cfg.zlim = [0,5];
   %cfg.zlim = [0 max(stati.stat,[],'all')];
   cfg.highlight        = 'on';
   int = all(stat_cmb.mask(:,idx(k,1):idx(k,2)), 2);
   % less severe
   % int = any(stat_cmb.mask(:,idx(k,1):idx(k,2)), 2); 
   cfg.highlightchannel = find(int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 12;
   cfg.colormap         = parula;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306cmb_helmet.mat';
   cfg.interactive      = 'no';
   cfg.parameter        = 'stat';
   ft_topoplotER(cfg, stat_cmb);
   set(gca,'ColorScale','linear') % if you want to plot on a logarithmic scale
   set(gca,'fontsize', 12)
end
posi        = get(subplot(3,5,13),'Position'); %[left bottom width height]
posi(1)     = posi(1)+posi(3)+0.01; posi(3) = 0.25*posi(3);
c           = colorbar('Position',posi);
c.LineWidth = 1;
c.FontSize  = 15;
title(c,'t-val','fontweight','bold');
%set(c,'Ticks',[0.1,1,2,3,4,5],'TickLabels',{'0.1','1','2','3','4','5'})

sgtitle('grandaverage - t-values - combined gradiometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_topoplot_timeseries_tvalues_combined_gradiometers.png'))

% timeseries
cfg               = [];
cfg.layout        = 'neuromag306cmb_helmet.mat';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor    = 'r';
figure('Position', get(0, 'Screensize'))
ft_multiplotER(cfg,stat_cmb);

sgtitle('grandaverage - timeseries combined gradiometers','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,'grandaverage_timeseries_combined_gradiometers.png'))



%% Clean up
addpath(fullfile('Z:','analysis','subject_files'));