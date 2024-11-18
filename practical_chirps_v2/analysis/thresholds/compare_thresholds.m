%--------------------------------------------------------------------------
% Till Habersetzer, 09.11.2024
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de
%--------------------------------------------------------------------------

close all
clearvars
clc

%% Load data
%--------------------------------------------------------------------------

current_dir = pwd;
cd('..')
settings_analysis
cd(current_dir)

% Load audiograms
%----------------
subjects     = [1,2,3];
n_sub        = length(subjects);
audiograms   = zeros(11,2,n_sub); % freq x left/right x subjects
subjectnames = cell(1,n_sub);

for sidx = 1:n_sub

    subidx = subjects(sidx);

    subject              = sprintf('sub-%02d', subidx);
    subjectnames{sidx}   = subject;
    data                 = readtable(fullfile(settings.path2project,'rawdata',subject,'beh', [subject '_task-audiogram_beh.tsv']), 'FileType', 'text', 'Delimiter', '\t');
    
    frequency            = data.frequency;
    audiograms(:,:,sidx) = [data.hearing_threshold_left,data.hearing_threshold_right];
    clear data

end

%% Load transient thresholds
%--------------------------------------------------------------------------
% Change from dB (p-p) peSPL to dB peSPL - correction factors
correction                              = struct();
correction.('tip300').left.click  = 5.5;
correction.('tip300').left.up     = 3.5;
correction.('tip300').right.click = 5;
correction.('tip300').right.up    = 3.5;

correction.('sensimetrics').left.click  = 3;
correction.('sensimetrics').left.up     = 3;
correction.('sensimetrics').right.click = 2.5;
correction.('sensimetrics').right.up    = 3.5;

apply_correction = false; % note: y-axis labels still in dB (p-p) peSPL


% Load Thresholds
%----------------
earphones    = {'tip300','sensimetrics'};
n_earphones  = length(earphones);
transients   = {'click','up'};
n_transients = length(transients);
sides        = {'left','right'};

thresholds              = struct();
thresholds.sensimetrics = zeros(n_sub,n_transients,2); % subjects x click/up x left/right
thresholds.tip300       = zeros(n_sub,n_transients,2);

for subidx= 1:n_sub % subject
    subject = subjectnames{subidx};

    for tidx = 1:n_transients % transient
        transient = transients{tidx};

        for sidx = 1:2 % side
            side = sides{sidx};

            for eidx = 1:n_earphones % earphone
                earphone = earphones{eidx};
    
                data = importdata(fullfile(settings.path2project,'rawdata',subject,'beh',sprintf('%s_stim-%s_%s_%s_hearinglevel_results.mat',subject,transient,side,earphone)));
                
                % Check
                if data.earphone ~= earphone
                    error('Unexpected earphone!')
                end
                if data.stimulus ~= transient
                    error('Unexpected transient!')
                end
                
                thresholds.(earphone)(subidx,tidx,sidx) = data.threshold_estimate;

                % Change from dB (p-p) peSPL to dB peSPL
                if apply_correction
                    thresholds.(earphone)(subidx,tidx,sidx) = thresholds.(earphone)(subidx,tidx,sidx) + correction.(earphone).(side).(transient);
                end
                clear data

            end

        end

    end

end

% Compute differences in threshold estimates
%-------------------------------------------
% Click - Up
idx_up    = find(contains(transients,'up'));
idx_click = find(contains(transients,'click'));

threshold_diff = struct();
for eidx = 1:n_earphones

    earphone  = earphones{eidx};
    threshold = thresholds.(earphone);

    thresholds_diff.(earphone) = squeeze(threshold(:,idx_click,:) - threshold(:,idx_up,:));

end

% Axis Limits
%------------
limits = zeros(2,2,2); % earphones x min/max x threshold+diff
for eidx = 1:n_earphones

    earphone  = earphones{eidx};
    limits(eidx,1,1) = min(thresholds.(earphone),[],'all');
    limits(eidx,2,1) = max(thresholds.(earphone),[],'all');

    limits(eidx,1,2) = min(thresholds_diff.(earphone),[],'all');
    limits(eidx,2,2) = max(thresholds_diff.(earphone),[],'all');
    
end

limits_min = squeeze(min(limits(:,1,:),[],1))-1; % minimum over both earphones 
limits_max = squeeze(max(limits(:,2,:),[],1))+1; 

% Visualize threshold
%--------------------
axs = cell(2,2);

for eidx = 1:n_earphones
    earphone       = earphones{eidx};
    threshold      = thresholds.(earphone);
    threshold_diff = thresholds_diff.(earphone);

    figure % new figure for each earphone
    sgtitle(earphone)

    for sidx = 1:2
        side = sides{sidx};
        if strcmp(side,'left')
            color = 'b';
        else
            color = 'r';
        end

        axs{1,sidx} = subplot(2,2,sidx); % threshold
        hold on
        axs{2,sidx} = subplot(2,2,2+sidx); % differences
        hold on

        for tidx = 1:n_transients

            transient = transients{tidx};
            if strcmp(transient,'up')
                marker = 'o';
            else
                marker = 'x';
            end
            % Thresholds
            p(sidx,tidx) = plot(axs{1,sidx},1:n_sub,threshold(:,tidx,sidx),'Color',color,'LineStyle','--','Marker',marker,'MarkerSize',15,'LineWidth',2);

            % Threshold differences
            bar(axs{2,sidx},1:n_sub,threshold_diff(:,sidx),'FaceColor','blue','BarWidth',0.25);

        end

    end
    set([axs{:}],'xtick', 1:n_sub, 'xticklabel',subjectnames,'xlim',[0.5,n_sub+0.5],'xgrid','on','ygrid','on')
    set([axs{1,:}],'ylim',[limits_min(1),limits_max(1)])
    set([axs{2,:}],'ylim',[limits_min(2),limits_max(2)])
    set([axs{1,1}],'Ylabel',ylabel('level / dB (p-p) peSPL'))
    set([axs{2,1}],'Ylabel',ylabel('level difference / dB'))
    xtickangle([axs{:}],45);
    title(axs{1,1},sides{1})
    title(axs{1,2},sides{2})
    title(axs{2,1},'Threshold: Click - Up')
    title(axs{2,2},'Threshold: Click - Up')
    legend(p(1,:),transients)

    clear threshold threshold_diff

end
    
%% Correlation with Pure-tone-average-4 (PTA4) 
%--------------------------------------------------------------------------
% average hearing thresholds for frequencies of 500, 1000, 2000, and 4000 Hz

idx  = ismember(frequency,[500,1000,2000,4000]);
pta4 = squeeze(mean(audiograms(idx,:,:),1))'; % subidx x (left/right)

% Plot
%-----
figure
plot(1:n_sub,pta4(:,1),'b--x','Markersize',15,'LineWidth',2)
hold on
plot(1:n_sub,pta4(:,2),'r--x','Markersize',15,'LineWidth',2)
ylabel('PTA4 / dB HL')
set(gca,'xtick', 1:n_sub, 'xticklabel',subjectnames)
xtickangle(45);
xlim([0.5,n_sub+0.5])
title('PTA4')
grid('on')
legend(sides)

for eidx = 1:n_earphones

    earphone = earphones{eidx};

    figure % new figure for each earphone
    sgtitle(earphone)

    axs = cell(1,2);
    for sidx = 1:2
        side = sides{sidx};
        if strcmp(side,'left')
            color = 'b';  
        else
            color = 'r';   
        end

        axs{sidx} = subplot(1,2,sidx);
        hold on

        for tidx = 1:n_transients

            transient = transients{tidx};
            if strcmp(transient,'up')
                marker = 'o';
            else
                marker = 'x';
            end

            p(tidx) = plot(pta4(:,sidx),thresholds.(earphone)(:,tidx,sidx),'Color',color,'MarkerSize',15,'LineStyle','none','Marker',marker,'LineWidth',2);
            
        end
        ylim([limits_min(1),limits_max(1)])
        xlim([min(pta4,[],'all')-0.5,max(pta4,[],'all')+0.5])
        grid('on')
        axis square
        title(side)
        xlabel('PTA4 / dB HL')
        ylabel('level / dB (p-p) peSPL')

    end
    legend(axs{1},transients)

end

%% Compare estimated thresholds and thresholds during the experiment
%--------------------------------------------------------------------------

% Load audiograms
%----------------
subjects     = [1,2,3];
n_sub        = length(subjects);
audiograms   = zeros(11,2,n_sub); % freq x left/right x subjects

thresholds_exp              = struct();
thresholds_exp.sensimetrics = zeros(n_sub,n_transients,2); % subjects x click/up x left/right
thresholds_exp.tip300       = zeros(n_sub,n_transients,2);
dB_thresholds               = struct;
for sidx = 1:n_sub
    subidx  = subjects(sidx);
    subject = sprintf('sub-%02d', subidx);

    for eidx = 1:n_earphones
        earphone = earphones{eidx};

        for tidx = 1:n_transients
            transient = transients{tidx};
            data      = importdata(fullfile(settings.path2project,'rawdata',subject,'meg', sprintf('%s_stim-%s_%s_results.mat',subject,transient,earphone)));
            
            if isfield(data,'target_level')
                thresholds_exp.(earphone)(subidx,tidx,:) = data.target_level - data.settings.dB_threshold;
            else
                thresholds_exp.(earphone)(subidx,tidx,:) = data.settings.calibration.(earphone).threshold;   
            end
            dB_thresholds.(earphone)(subidx,tidx,:) = data.settings.dB_threshold; % just store for later, dB SL computation
            
            clear data
        end
    end
end

% Compute differences in threshold estimates and thresholds during
% experiment
%-----------------------------------------------------------------
% Click - Up

thresholds_exp_diff = struct();
levels_exp          = struct(); % levels in dB SL during experiment
for eidx = 1:n_earphones

    earphone     = earphones{eidx};
    dB_threshold = dB_thresholds.(earphone);

    thresholds_exp_diff.(earphone) = thresholds_exp.(earphone)-thresholds.(earphone);
    levels_exp.(earphone)          = thresholds_exp.(earphone)-thresholds.(earphone) + dB_threshold; 

end

% Axis Limits
%------------
limits = zeros(2,2,3); % earphones x min/max x threshold+diff+leveldBSL
for eidx = 1:n_earphones

    earphone         = earphones{eidx};
    limits(eidx,1,1) = min([thresholds.(earphone),thresholds_exp.(earphone)],[],'all');
    limits(eidx,2,1) = max([thresholds.(earphone),thresholds_exp.(earphone)],[],'all');

    limits(eidx,1,2) = min(thresholds_exp_diff.(earphone),[],'all');
    limits(eidx,2,2) = max(thresholds_exp_diff.(earphone),[],'all');

    limits(eidx,1,3) = min(levels_exp.(earphone),[],'all');
    limits(eidx,2,3) = max(levels_exp.(earphone),[],'all');
    
end

limits_min = squeeze(min(limits(:,1,:),[],1))-1; % minimum over both earphones 
limits_max = squeeze(max(limits(:,2,:),[],1))+1; 

% Visualize threshold
%--------------------
axs          = cell(3,2);
legend_names = cell(n_transients,2);

for eidx = 1:n_earphones
    earphone           = earphones{eidx};
    threshold          = thresholds.(earphone);
    threshold_exp      = thresholds_exp.(earphone);
    threshold_exp_diff = thresholds_exp_diff.(earphone);
    level_exp          = levels_exp.(earphone);

    figure % new figure for each earphone
    sgtitle(earphone)

    for sidx = 1:2
        side = sides{sidx};
        if strcmp(side,'left')
            color1 = 'b';
            color2 = [0.3010 0.7450 0.9330];
        else
            color1 = 'r';
            color2 = [0.8500 0.3250 0.0980];
        end

        axs{1,sidx} = subplot(3,2,sidx); % threshold
        hold on
        axs{2,sidx} = subplot(3,2,2+sidx); % differences transients
        hold on
        axs{3,sidx} = subplot(3,2,4+sidx); % level in dB SL
        hold on

        for tidx = 1:n_transients

            transient = transients{tidx};
            if strcmp(transient,'up')
                marker = 'o';
            else
                marker = 'x';
            end
            % Thresholds
            p1(sidx,tidx) = plot(axs{1,sidx},1:n_sub,threshold(:,tidx,sidx),'Color',color1,'LineStyle','--','Marker',marker,'MarkerSize',15,'LineWidth',2);
            p2(sidx,tidx) = plot(axs{1,sidx},1:n_sub,threshold_exp(:,tidx,sidx),'Color',color2,'LineStyle','--','Marker',marker,'MarkerSize',15,'LineWidth',2);

            % Threshold differences
            b = bar(axs{2,sidx},1:n_sub,threshold_exp_diff(:,:,sidx),'grouped','FaceColor','blue','BarWidth',0.25);
            b(1).FaceColor = [0.2 0.2 0.8]; 
            b(2).FaceColor = [0.8 0.2 0.2]; 

            legend_names{tidx,1} = sprintf("%s (threshold)",transient);
            legend_names{tidx,2} = sprintf("%s (during experiment)",transient);

            % Level in dB SL
            c              = bar(axs{3,sidx},1:n_sub,level_exp(:,:,sidx),'grouped','FaceColor','blue','BarWidth',0.25);
            c(1).FaceColor = [0.2 0.2 0.8]; 
            c(2).FaceColor = [0.8 0.2 0.2]; 

            data = level_exp(:,:,sidx);
            for i = 1:size(data, 2)
                % Get the x-coordinates of the bars for this condition
                x = b(i).XEndPoints;
                % Get the y-values (heights) of the bars
                y = data(:, i);
                % Add the text on top of the bars
                text(x, y, string(y), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10);
            end
            clear data

        end
          
    end

    set([axs{:}],'xtick', 1:n_sub, 'xticklabel',subjectnames,'xlim',[0.5,n_sub + 1],'xgrid','on','ygrid','on')
    set([axs{1,:}],'ylim',[limits_min(1),limits_max(1)])
    set([axs{2,:}],'ylim',[limits_min(2),limits_max(2)])
    set([axs{3,:}],'ylim',[0,limits_max(3)+5])
    set([axs{1,1}],'Ylabel',ylabel('level / dB (p-p) peSPL'))
    set([axs{2,1}],'Ylabel',ylabel('level difference / dB'))
    set([axs{3,1}],'Ylabel',ylabel('level / dB SL'))
    xtickangle([axs{:}],45);
    title(axs{1,1},sides{1})
    title(axs{1,2},sides{2})
    title(axs{2,1},'Threshold differences (experiment - measured)')
    title(axs{2,2},'Threshold differenses (experiment - measured)')
    title(axs{3,1},'Level during experiment in dB SL')
    title(axs{3,2},'Level during experiment in dB SL')
    legend(axs{1,1},[p1(1,:),p2(1,:)],[legend_names{:,1},legend_names{:,2}])
    legend(axs{1,2},[p1(2,:),p2(2,:)],[legend_names{:,1},legend_names{:,2}])
    legend(axs{2,1},transients)
    legend(axs{3,1},transients)

    clear threshold threshold_diff

end



