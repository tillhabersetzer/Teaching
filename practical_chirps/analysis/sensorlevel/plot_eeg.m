close all
clear 
clc 

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

conditions = settings.conditions;

%% Compute N19P30 32 channel EEG CAP
%--------------------------------------------------------------------------
 
subidx  = 1;
subject = ['sub-',num2str(subidx,'%02d')];
data    = importdata(fullfile(settings.path2project,'derivatives',subject,'sensorlevel',[subject,'_erp-N19P30.mat']));

avg       = data.avgs;
layout    = data.layout;
reference = data.reference;
clear data

cfg        = [];
cfg.layout = layout;
ft_layoutplot(cfg)

% topoplot
%---------
cfg            = [];
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbk';
cfg.showlabels = 'yes';
cfg.comment    = 'nose reference';
ft_multiplotER(cfg, avg{1,1}, avg{1,2}, avg{1,3});
legend({conditions{1};conditions{2};conditions{3}});
title('nose reference')

cfg            = [];
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbk';
cfg.showlabels = 'yes';
cfg.comment    = 'avg reference';
ft_multiplotER(cfg, avg{2,1}, avg{2,2}, avg{2,3});
legend({conditions{1};conditions{2};conditions{3}});
title('avg reference')

% single channel
%---------------
chan2plot = 'EEG005';
cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg{2,:});
legend({conditions{1};conditions{2};conditions{3}});
title([subject,': ', chan2plot,' (avg reference)'])

%% Visualize N19P30 10 channel cEEGrid
%--------------------------------------------------------------------------

subidx  = 2;
subject = ['sub-',num2str(subidx,'%02d')];
data    = importdata(fullfile(settings.path2project,'derivatives',subject,'sensorlevel',[subject,'_erp-N19P30_cEEGrid.mat']));

avgs_left          = data.avgs_left;
avgs_right         = data.avgs_right;
avgs_between       = data.avgs_between;
layout             = data.layout;
combinations_label = data.combinations_label;
clear data

cfg        = [];
cfg.layout = layout;
ft_layoutplot(cfg)

% limits
mini = 1;
maxi = 0;
for c=1:10
    for cidx=1:3
        a = min([avgs_left{c,cidx}.avg,avgs_right{c,cidx}.avg]);
        b = max([avgs_left{c,cidx}.avg,avgs_right{c,cidx}.avg]);
        if a<mini
            mini = a;
        end
        if b>maxi
            maxi = b;
        end
    end
end

% add channel combinations between ears
for c=1:25
    for cidx=1:3
        a = min([avgs_between{c,cidx}.avg]);
        b = max([avgs_between{c,cidx}.avg]);
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

% right ear
%----------
figure
for c=1:10 % loop over channel combinations
    subplot(2,5,c);
    for cidx = 1:3 % loop over conditions
        hold on
        plot(avgs_right{c,cidx}.time*1000,avgs_right{c,cidx}.avg,color{cidx});
    end
    xlabel('t/ms')
    
    if c==1
        legend({conditions{1};conditions{2};conditions{3}});
    end
    axis(axisvec) 
    title([combinations_label{1}{c,1},'-',combinations_label{1}{c,2}])
end
sgtitle([subject,': right ear'])

% left ear
%----------
figure
for c=1:10 % loop over channel combinations
    subplot(2,5,c);
    for cidx = 1:3 % loop over conditions
        hold on
        plot(avgs_left{c,cidx}.time*1000,avgs_left{c,cidx}.avg,color{cidx});
    end
    xlabel('t/ms')
    
    if c==1
        legend({conditions{1};conditions{2};conditions{3}});
    end
    axis(axisvec) 
    title([combinations_label{2}{c,1},'-',combinations_label{2}{c,2}])
end
sgtitle([subject,': left ear'])

% between ears
%-------------
figure
for c=1:25 % loop over channel combinations
    subplot(5,5,c);
    for cidx = 1:3 % loop over conditions
        hold on
        plot(avgs_between{c,cidx}.time*1000,avgs_between{c,cidx}.avg,color{cidx});
    end
    xlabel('t/ms')
    
    if c==1
        legend({conditions{1};conditions{2};conditions{3}});
    end
    axis(axisvec) 
    title([combinations_label{3}{c,1},'-',combinations_label{3}{c,2}])
end
sgtitle([subject,': between ears'])

%% Missmatch negativity 10 channel cEEGrid
%--------------------------------------------------------------------------

subidx  = 2;
subject = ['sub-',num2str(subidx,'%02d')];
data    = importdata(fullfile(settings.path2project,'derivatives',subject,'sensorlevel',[subject,'_erp-mmn_cEEGrid.mat']));

avgs_left          = data.avgs_left;
avgs_right         = data.avgs_right;
avgs_between       = data.avgs_between;
layout             = data.layout;
combinations_label = data.combinations_label;
order              = data.order;
clear data

cfg        = [];
cfg.layout = layout;
ft_layoutplot(cfg)

% limits
mini = 1;
maxi = 0;
for c=1:10
    for cidx=1:3
        a = min([avgs_left{c,cidx}.avg,avgs_right{c,cidx}.avg]);
        b = max([avgs_left{c,cidx}.avg,avgs_right{c,cidx}.avg]);
        if a<mini
            mini = a;
        end
        if b>maxi
            maxi = b;
        end
    end
end

% add channel combinations between ears
for c=1:25
    for cidx=1:3
        a = min([avgs_between{c,cidx}.avg]);
        b = max([avgs_between{c,cidx}.avg]);
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
color   = {'r','b','k'};

% right ear
%----------
figure
for c=1:10 % loop over channel combinations
    subplot(2,5,c);
    for cidx = 1:3 % loop over conditions
        hold on
        plot(avgs_right{c,cidx}.time*1000,avgs_right{c,cidx}.avg,color{cidx});
    end
    xlabel('t/ms')
    
    if c==1
        legend({order{1};order{2};order{3}});
    end
    axis(axisvec) 
    title([combinations_label{1}{c,1},'-',combinations_label{1}{c,2}])
end
sgtitle([subject,': right ear'])

% left ear
%----------
figure
for c=1:10 % loop over channel combinations
    subplot(2,5,c);
    for cidx = 1:3 % loop over conditions
        hold on
        plot(avgs_left{c,cidx}.time*1000,avgs_left{c,cidx}.avg,color{cidx});
    end
    xlabel('t/ms')
    
    if c==1
        legend({order{1};order{2};order{3}});
    end
    axis(axisvec) 
    title([combinations_label{2}{c,1},'-',combinations_label{2}{c,2}])
end
sgtitle([subject,': left ear'])

% between ears
%-------------
figure
for c=1:25 % loop over channel combinations
    subplot(5,5,c);
    for cidx = 1:3 % loop over conditions
        hold on
        plot(avgs_between{c,cidx}.time*1000,avgs_between{c,cidx}.avg,color{cidx});
    end
    xlabel('t/ms')
    
    if c==1
        legend({order{1};order{2};order{3}});
    end
    axis(axisvec) 
    title([combinations_label{3}{c,1},'-',combinations_label{3}{c,2}])
end
sgtitle([subject,': between ears'])

%% plot cEEGrid sensor positions

subidx   = 2;
subject  = ['sub-',num2str(subidx,'%02d')];

datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-clicks.fif']);
sens      = ft_read_sens(datapath,'senstype','eeg');
sens      = ft_convert_units(sens,'mm');
mri       = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_T1w-segmented.mat']));
headshape = ft_read_headshape(datapath,'unit','mm');

% scalp surface
cfg             = [];
cfg.tissue      = {'scalp'};
cfg.method      = 'isosurface';
cfg.numvertices = 10000;
scalp           = ft_prepare_mesh(cfg, mri);

% manipulate sens for sub-01
if strcmp(subject,'sub-01')
    sens.chanpos  = sens.chanpos(33:42,:);
    sens.chantype = sens.chantype(33:42);
    sens.chanunit = sens.chanunit(33:42);
    sens.elecpos  = sens.elecpos(33:42,:);
    sens.label    = sens.label(33:42);

    for c = 1:10
        sens.label{c} = ['EEG0',num2str(c,'%02d')];
    end
end

figure
ft_plot_mesh(headshape.fid, 'vertexcolor', 'g', 'vertexsize', 10);
ft_plot_mesh(scalp, 'facecolor', 'skin')
%ft_plot_headshape(headshape, 'vertexcolor', 'k');
ft_plot_sens(sens,'label' ,'on','coilshape','disc','facecolor','r')
lighting phong
camlight left
camlight right
material dull
alpha 0.5
sgtitle(subject)


%% Clean-Up
rmpath('..')

