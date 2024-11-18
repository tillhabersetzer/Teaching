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

% apply ica
ica_status = 1;

check = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'));

if ica_status
    add = '_ica';
else 
    add = '';
end

% load preprocessed data
%-----------------------
N_subj  = length(subjects);
avg_pre = cell(1,N_subj);
avg_pst = cell(1,N_subj);
info    = cell(1,N_subj);
trials  = zeros(1,N_subj);

for s = 1:N_subj
    eval(subjects{s})
    data            = importdata(fullfile(subjectdata.chirps_sensorlevel,[subjectdata.subjectname,'_averages',add,'.mat']));
    avg             = data.avg_w;
    cfg             = [];
    cfg.latency     = [-0.45 0];
    avg_pre{s}      = ft_selectdata(cfg,avg);
    cfg.latency     = [0 0.45];
    avg_pst{s}      = ft_selectdata(cfg,avg);
    avg_pre{s}.time = avg_pst{s}.time; % change time axis
    clear avg
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

%% calculate grandaverage 
%------------------------
cfg         = [];
cfg.latency = 'all';
gavg_pre    = ft_timelockgrandaverage(cfg,avg_pre{:});
gavg_pst    = ft_timelockgrandaverage(cfg,avg_pst{:});

cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffect     = ft_math(cfg,gavg_pst,gavg_pre);

if check
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg,gavg_pre,gavg_pst);
    sgtitle('grandaverage: all channels') 
    
    figure
    ft_multiplotER(cfg,raweffect);
    sgtitle('raweffect: all channels')   
    
    % check individual averages
    figure
    ft_multiplotER(cfg,avg_pre{10},avg_pst{10});
    sgtitle('raweffect: all channels')   
end

%% statistic
%-----------
cfg                  = [];
cfg.channel          = 'meg';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 4;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
cfg_neighb           = [];
cfg_neighb.channel   = 'meg';
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, avg_pre{1}.grad);
%ft_neighbourplot(cfg,datapre_w)

design                      = zeros(2,2*N_subj);
design(1,1:N_subj)          = 1;
design(1,N_subj+1:2*N_subj) = 2;
design(2,1:N_subj)          = [1:N_subj];
design(2,N_subj+1:2*N_subj) = [1:N_subj];

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable 
cfg.uvar   = 2; % row of design matrix that contains unit of observation

stat = ft_timelockstatistics(cfg,avg_pst{:},avg_pre{:});

%% make nice plots
%-----------------
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust  = find(pos_cluster_pvals < stat.cfg.alpha); % use all significant clusters
pos               = ismember(stat.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < stat.cfg.alpha);
neg               = ismember(stat.negclusterslabelmat, neg_signif_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffect.label, stat.label);

% time windows for topoplots
windows = [0,0.05;0.05,0.075;0.075,0.1;0.1,0.125;0.125,0.15;0.15,0.175;0.175,0.2;...
           0.2,0.225;0.225,0.25;0.25,0.3;0.3,0.35;0.35,0.4;0.4,0.45];
idx     = dsearchn(stat.time',reshape(windows,1,numel(windows))');
idx     = reshape(idx,size(windows,1),size(windows,2));

% if check
    
new_dir = fullfile('Z:','analysis','generated_data','grouplevel','chirps','sensorlevel');
if ~exist(new_dir, 'dir')
   mkdir(new_dir)
end

% Plot over raweffect
%--------------------
figure('Position', get(0, 'Screensize'))
for k = 1:13
   subplot(3,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   pos_int  = zeros(numel(raweffect.label),1);
   neg_int  = zeros(numel(raweffect.label),1);
   %pos_int(i1) = all(pos(i2,idx(k,1):idx(k,2)), 2);
   %neg_int(i1) = all(neg(i2,idx(k,1):idx(k,2)), 2);
   % less severe
   pos_int(i1) = any(pos(i2,idx(k,1):idx(k,2)), 2);
   neg_int(i1) = any(neg(i2,idx(k,1):idx(k,2)), 2);
   cfg.highlight        = 'on';
   % Get the index of each significant channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 20;
   cfg.highlightsymbol  = '*';
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306all_helmet.mat';
   cfg.interactive      = 'no';
   ft_topoplotER(cfg, raweffect);
   set(gca,'fontsize', 12)
end
sgtitle('grandaverage - raweffect','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,['grandaverage_topoplot_timeseries_raweffect',add,'.png']))

% Plot over t-values
%-------------------
% change t-values to abs for plotting
stati      = stat;
stati.stat = abs(stat.stat);

figure('Position', get(0, 'Screensize'))
for k = 1:10
   subplot(2,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   cfg.zlim = [0,5];
   %cfg.zlim = [0 max(stati.stat,[],'all')];
   cfg.highlight        = 'on';
   int = all(stati.mask(:,idx(k,1):idx(k,2)), 2);
   % less severe
%    int = any(stati.mask(:,idx(k,1):idx(k,2)), 2); 
   cfg.highlightchannel = find(int);
   cfg.highlightcolor   = [1 0 0];
   cfg.highlightsize    = 20;
   cfg.highlightsymbol  = '*';
   cfg.colormap         = parula;
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = 'neuromag306all_helmet.mat';
   cfg.interactive      = 'no';
   cfg.parameter        = 'stat';
   ft_topoplotER(cfg, stati);
   set(gca,'ColorScale','log') % if you want to plot on a logarithmic scale
   set(gca,'fontsize', 12)
end
posi        = get(subplot(2,5,10),'Position'); %[left bottom width height]
posi(1) = posi(1)+posi(3)+0.01; posi(3) = 0.25*posi(3);
c           = colorbar('Position',posi);
c.LineWidth = 1;
c.FontSize  = 15;
title(c,'|t-val|','fontweight','bold');
set(c,'Ticks',[0.1,1,2,3,4,5],'TickLabels',{'0.1','1','2','3','4','5'})

sgtitle('grandaverage - t-values','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,['grandaverage_topoplot_timeseries_tvalues',add,'.png']))

% timeseries
cfg               = [];
cfg.layout        = 'neuromag306all_helmet.mat';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor    = 'r';
figure('Position', get(0, 'Screensize'))
ft_multiplotER(cfg,stat);

sgtitle('grandaverage - timeseries','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,['grandaverage_timeseries',add,'.png']))

% end

%% Clean up
addpath(fullfile('Z:','analysis','subject_files'));