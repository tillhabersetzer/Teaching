close all; clear all; clc;

% Settings
%--------------------------------------------------------------------------
% choose subject 
% subjectlist = {'subject02'};
for i = 2:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjectlist{i-1} = [subject,num2str(i)]; 
end

% choose files
files2preproc = 'stories_maxfilter';

% choose envelope type for crosscorrelation
% envelopetype = 'onset_envelope';
envelopetype = 'envelope';

% load data with performed ica (1) or without (0)
ica_on = 1; 

% load data with whitening (1) or without (0)
% not necessary for correlation in sensorspace! cross correlation is
% normalized
whitening = 0; 

% load data with additional zscoring of all trials
zscoring = 1;

% check data
check = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'))

if ica_on
    add = '_ica';
else 
    add = '';
end

if zscoring
    add2 = '_zscored';
else 
    add2 = '';
end

new_dir = fullfile('Z:','analysis','generated_data','grouplevel','speech',...
          'crosscorrelation',envelopetype,'sensorlevel');
if ~exist(new_dir, 'dir')
   mkdir(new_dir)
end

%% load data
%--------------------------------------------------------------------------
N_subj      = length(subjectlist);
avg         = cell(1,N_subj);
avg_shuffle = cell(1,N_subj);
info        = cell(1,N_subj);

for s = 1:N_subj
    eval(subjectlist{s})
    if ica_on
       filename_new = [subjectdata.subjectname '_' files2preproc '_crosscorr_ica'];
    else
       filename_new = [subjectdata.subjectname '_' files2preproc '_crosscorr_noica'];
    end
    if whitening
        filename_new = [filename_new '_whitened'];
    end
    if zscoring
        filename_new = [filename_new '_zscored'];
    end
    % data
    %-----
    data           = importdata(fullfile(subjectdata.speech,'crosscorrelation',...
                     envelopetype,'sensorlevel',[filename_new,'.mat']));
    avg{s}         = data.avg;
    avg_shuffle{s} = data.avg_shuffle;
    info{s}        = data.info;  
    clear data
    disp([subjectlist{s},' loaded.'])
end

%% grandaverage
%--------------------------------------------------------------------------
cfg              = [];
cfg.latency      = 'all';
grandavg         = ft_timelockgrandaverage(cfg,avg{:});
grandavg_shuffle = ft_timelockgrandaverage(cfg,avg_shuffle{:});

cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffect     = ft_math(cfg,grandavg,grandavg_shuffle);

if check
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg,grandavg,grandavg_shuffle);
    sgtitle('grandaverage: all channels') 
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg,raweffect);
    sgtitle('grandaverage raweffect: all channels') 
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306mag_helmet.mat';
    ft_multiplotER(cfg,grandavg,grandavg_shuffle);
    sgtitle('grandaverage: magnetometer') 
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306planar_helmet.mat';
    ft_multiplotER(cfg,grandavg,grandavg_shuffle);
    sgtitle('grandaverage: gradiometer') 
    
    % choose subject to look at
    n              = 18;
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306all_helmet.mat';
    ft_multiplotER(cfg,avg{n},avg_shuffle{n});
    sgtitle([subjectlist{n},' | average: gradiometer']) 
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
cfg.minnbchan        = 5;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb           = [];
cfg_neighb.channel   = 'meg';
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, avg{1}.grad);
%ft_neighbourplot(cfg,datapre_w)

design                      = zeros(2,2*N_subj);
design(1,1:N_subj)          = 1;
design(1,N_subj+1:2*N_subj) = 2;
design(2,1:N_subj)          = [1:N_subj];
design(2,N_subj+1:2*N_subj) = [1:N_subj];

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable 
cfg.uvar   = 2; % row of design matrix that contains unit of observation

stat = ft_timelockstatistics(cfg,avg{:},avg_shuffle{:});

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
[i1,i2] = match_str(raweffect.label,stat.label);

% time windows for topoplots
windows = [-0.1,0;0,0.1;0.1,0.2;0.2,0.3;0.3,0.4;0.4,0.5;0.5,0.6;...
           0.6,0.7;0.7,0.8;0.8,0.9];
idx     = dsearchn(stat.time',reshape(windows,1,numel(windows))');
idx     = reshape(idx,size(windows,1),size(windows,2));

if check
   
% Plot over raweffect
%--------------------
figure('Position', get(0, 'Screensize'))
for k = 1:10
   subplot(2,5,k);
   cfg      = [];
   cfg.xlim = windows(k,:);  
   pos_int  = zeros(numel(raweffect.label),1);
   neg_int  = zeros(numel(raweffect.label),1);
   pos_int(i1) = all(pos(i2,idx(k,1):idx(k,2)), 2);
   neg_int(i1) = all(neg(i2,idx(k,1):idx(k,2)), 2);
   % less severe
%    pos_int(i1) = any(pos(i2,idx(k,1):idx(k,2)), 2);
%    neg_int(i1) = any(neg(i2,idx(k,1):idx(k,2)), 2);
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
   ft_topoplotER(cfg,raweffect);
   set(gca,'fontsize', 12)
end
sgtitle('grandaverage - raweffect','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,['grandaverage_topoplot_crosscorrelation_raweffect',add,add2,'.png']))

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
   cfg.zlim = [0,8];
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
set(c,'Ticks',[0.1,1,2,3,4,5,8],'TickLabels',{'0.1','1','2','3','4','5','8'})

sgtitle('grandaverage - t-values','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,['grandaverage_topoplot_crosscorrelation_tvalues',add,add2,'.png']))

% timeseries
%-----------
cfg               = [];
cfg.layout        = 'neuromag306all_helmet.mat';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor    = 'r';
figure('Position', get(0, 'Screensize'))
ft_multiplotER(cfg,stat);

sgtitle('grandaverage - crosscorrelation','fontweight','bold','FontSize',20)
saveas(gcf,fullfile(new_dir,['grandaverage_crosscorrelation',add,add2,'.png']))


end


%% Clean up
rmpath(fullfile('Z:','analysis','subject_files'))
