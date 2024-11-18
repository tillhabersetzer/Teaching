close all
clear 
clc 

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','..'))
eval('main_settings')

% check whether maxfiltered data should be analyzed
maxfilter = settings.maxfilter;
% acoustic latency correction for triggers
latency_correction = settings.latency_correction;

%% Script settings
%--------------------------------------------------------------------------
subidx  = 1;
%--------------------------------------------------------------------------

% conditions
conditions = settings.conditions;

% Trigger IDs
TrigID.click     = 1; 
TrigID.upchirp   = 2; 
TrigID.downchirp = 4;

%% Preprocess data
%--------------------------------------------------------------------------
subject           = ['sub-',num2str(subidx,'%02d')];
data_preprocessed = cell(1,3);
avg               = cell(1,3);
N_trials          = zeros(1,3);

for cidx=1:3 % loop over conditions

    % Path for MEG data
    %------------------------------------------------------------------
    condition = conditions{cidx};

    if maxfilter
        datapath = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-',condition,'-raw_tsss.fif']);
    else
        datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);
    end

    switch conditions{cidx}
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
    cfg.trialfun            = 'ft_trialfun_general'; 
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
    cfg.channel      = 'meg'; 
    cfg.continuous   = 'yes';
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = [16,120];
    cfg.dftfilter    = 'yes';       % enable notch filtering to eliminate power line noise
    cfg.dftfreq      = [50,100];    % set up the frequencies for notch filtering
    cfg.dftreplace   = 'zero';      % 'zero' is default'
    cfg.coilaccuracy = 0;           % ensure that sensors are expressed in SI units
    data             = ft_preprocessing(cfg); 

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

    data_preprocessed{cidx} = epoched;

    % Note down amount of preserved epochs
    %-------------------------------------
    N_trials(cidx) = length(epoched.trial);

    % Timelockanalysis
    %-----------------
    cfg                  = [];
    cfg.keeptrials       = 'yes';
    cfg.covariancewindow = 'all';
    cfg.covariance       = 'yes';
    avg{cidx}            = ft_timelockanalysis(cfg, epoched);

    clear data epoched

end % loop over conditions

% Compute pooled average (over clicks and downchirps)
%----------------------------------------------------
idx              = find(contains(conditions,{'clicks','downchirps'}));

cfg              = [];
data_pooled      = ft_appenddata(cfg, data_preprocessed{idx(1)}, data_preprocessed{idx(2)});
data_pooled.grad = data_preprocessed{1}.grad;

cfg        = [];
avg_pooled = ft_timelockanalysis(cfg, data_pooled);

%% Add latency correction 
%--------------------------------------------------------------------------
% ~5ms 
% 1.8ms transmission delay of sound + max 3ms Recording of analog MEG signals
% -> digital triggers arrive 5 ms to early in MEG recordings

for cidx=1:3
    avg{cidx}.time = avg{cidx}.time - latency_correction; % in sec
end
avg_pooled.time = avg_pooled.time - latency_correction;

%% Visualize results
layout = 'neuromag306mag.lay';
% layout = 'neuromag306planar.lay';

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
cfg.linecolor  = 'rbk';
ft_multiplotER(cfg, avg{:});
legend({conditions{1};conditions{2};conditions{3}});
title(subject)
%% Source Analysis

% Switch between erf-types
%-------------------------
timewin  = [0.02 0.05]; % 20-50ms covering N19m-P30m-Peak2Peak interval
chantype = 'megplanar'; % 'megmag','megplanar','meg'

% Load sourcemodel and headmodel
%-------------------------------
% Proceed in SI-units. Leads apparantly to correct units for dipole moments in Am

sourcemodel = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_sourcemodel-volumetric.mat']));
sourcemodel = ft_convert_units(sourcemodel, 'm'); 

headmodel = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_headmodel-singleshell.mat']));
headmodel = ft_convert_units(headmodel, 'm');  


%% Fit initial dipole model to pooled conditions (clicks+downchirps)
%----------------------------------------------------------------------       

% Apply Sphering if magnetometers and gradiometers are used
% together
%----------------------------------------------------------
sphering = strcmp(chantype,'meg');
if sphering     
    if maxfilter
        noisefile = fullfile(settings.path2project,'derivatives',subject,'maxfilter',[subject,'_task-emptyroom-raw_tsss.fif']);
    else
        noisefile = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-emptyroom.fif']);
    end

    cfg              = [];
    cfg.dataset      = noisefile;
    cfg.channel      = 'meg';
    cfg.continuous   = 'yes';
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = [16,120];
    cfg.dftfilter    = 'yes';     
    cfg.dftfreq      = [50,100];    
    cfg.dftreplace   = 'zero';
    cfg.coilaccuracy = 0; 
    noise            = ft_preprocessing(cfg);

    cfg                  = [];
    cfg.removemean       = 'yes'; % default for covariance computation
    cfg.covariance       = 'yes';
    cfg.covariancewindow = 'all';
    avg_noise            = ft_timelockanalysis(cfg,noise);

%     imagesc(avg_noise.cov)
end

% Symmetrical gridsearch
%-----------------------
cfg             = [];
cfg.latency     = timewin;
cfg.numdipoles  = 2; % fit 2 dipoles
cfg.symmetry    = 'x'; % expect activity in both auditory cortices, fit symmetric dipole
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel;
cfg.gridsearch  = 'yes';
cfg.senstype    = 'meg';
cfg.channel     = chantype;

if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
source_pooled_sym = ft_dipolefitting(cfg, avg_pooled);

% Now that we have a better starting point for the dipole fit, we can release the symmetry contstraint
%-----------------------------------------------------------------------------------------------------
cfg                  = [];
cfg.latency          = timewin;
cfg.numdipoles       = 2;
cfg.symmetry         = []; % no symmetry constraint
cfg.gridsearch       = 'no';
cfg.nonlinear        = 'yes'; % non-linear search for optimal dipole parameters
cfg.dip.pos          = source_pooled_sym.dip.pos;
cfg.sourcemodel.unit = sourcemodel.unit; % defines units of dipole positions
cfg.headmodel        = headmodel;
cfg.senstype         = 'meg';
cfg.channel          = chantype;

if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
source_pooled_nosym = ft_dipolefitting(cfg, avg_pooled);

%% Beamforming

% Computation of leadfield for predefined dipole positions
%--------------------------------------------------------------------------
cfg                       = [];
cfg.grad                  = avg_pooled.grad;  % gradiometer distances
cfg.channel               = chantype;                       
cfg.sourcemodel.pos       = source_pooled_nosym.dip.pos;
cfg.sourcemodel.unit      = sourcemodel.unit;
cfg.headmodel             = headmodel;                    
lf                        = ft_prepare_leadfield(cfg);

% Common filter
%--------------------------------------------------------------------------
cfg      = [];
data_all = ft_appenddata(cfg, data_preprocessed{:});
data_all.grad = data_preprocessed{1}.grad;

cfg                  = [];
cfg.covariancewindow = 'all';
cfg.covariance       = 'yes';
avg_all              = ft_timelockanalysis(cfg, data_all);

% create spatial filter using the lcmv beamformer
% check out for more details: ft_inverse_lcmv
cfg                 = [];
cfg.channel         = chantype;
cfg.method          = 'lcmv';
cfg.sourcemodel     = lf; % leadfield
cfg.headmodel       = headmodel; % volume conduction model (headmodel)
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes'; % project on axis of most variance using SVD
source_all          = ft_sourceanalysis(cfg, avg_all);

% project all trials through common spatial filter 
source = cell(1,3);

cfg                    = [];
cfg.channel            = chantype;
cfg.method             = 'lcmv';
cfg.sourcemodel        = lf;      
cfg.headmodel          = headmodel;       
cfg.sourcemodel.filter = source_all.avg.filter; 
source{1}              = ft_sourceanalysis(cfg, avg{1}); 
source{2}              = ft_sourceanalysis(cfg, avg{2}); 
source{3}              = ft_sourceanalysis(cfg, avg{3}); 

%% Visualize results

% Mapping between dipole locations and hemispheres
%-------------------------------------------------
pos     = source_all.pos;
mapping = check_diploc(pos); 

name = {'left hemisphere','right hemisphere'};
time = source_all.time*1000; 

mini = min([source{1}.avg.mom{1},source{1}.avg.mom{2},source{2}.avg.mom{1},source{2}.avg.mom{1} ...
            source{3}.avg.mom{1},source{3}.avg.mom{2}],[],'all');
maxi = max([source{1}.avg.mom{1},source{1}.avg.mom{2},source{2}.avg.mom{1},source{2}.avg.mom{1} ...
            source{3}.avg.mom{1},source{3}.avg.mom{2}],[],'all');
time2plot = [-50 300];
axisvec   = horzcat(time2plot,[mini,maxi]);

figure
for i=1:2
    idx  = mapping(i);
       
    subplot(2,1,i); 
    hold on
    plot(time, source{1}.avg.mom{idx}, '-')
    plot(time, source{2}.avg.mom{idx}, '-')
    plot(time, source{3}.avg.mom{idx}, '-')
    legend({'click', 'upchirp', 'downchirp'});
    axis(axisvec)
    grid on
    title(name{i})
    if i ==2
        xlabel('t / ms')
    end
%     ylabel('dipole moment / nAm')
end
sgtitle([subject,': fixed oriented dipole via beamformer'])


%% Beamforming single condition

% Computation of leadfield for predefined dipole positions
%--------------------------------------------------------------------------
chantype = 'megmag';

cfg                       = [];
cfg.grad                  = avg{1}.grad;  % gradiometer distances
cfg.channel               = chantype;                       
cfg.sourcemodel           = sourcemodel;
cfg.headmodel             = headmodel;     
cfg.singleshell.batchsize = 2000;
% cfg.normalize             = 'yes'; 
lf                        = ft_prepare_leadfield(cfg);

% create spatial filter using the lcmv beamformer
% check out for more details: ft_inverse_lcmv
cfg                 = [];
cfg.channel         = chantype;
cfg.method          = 'lcmv';
cfg.sourcemodel     = lf; % leadfield
cfg.headmodel       = headmodel; % volume conduction model (headmodel)
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.kappa      = 70;
% cfg.lcmv.lambda     = '5%';
% cfg.lcmv.weightnorm = 'unitnoisegain';
% cfg.lcmv.weightnorm = 'nai';
cfg.lcmv.projectnoise = 'yes';
source1             = ft_sourceanalysis(cfg, avg{1});

% NAI
source1.avg.pow = source1.avg.pow ./ source1.avg.noise;

% interpolate the source power onto the individual MRI
mri = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_T1w-segmented.mat']));
% Keep it in SI-units
mri = ft_convert_units(mri, 'm'); 
cfg              = [];
cfg.parameter    = 'avg.pow';
cfg.voxelcoord   = 'no';
cfg.interpmethod = 'nearest';
source1_int      = ft_sourceinterpolate(cfg, source1, mri);

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
% cfg.location = [-42 -18 67];
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, source1_int);

%% Condition contrast
chantype = 'megmag';

cfg        = [];
cfg.toilim = [-0.05 0]; % does not make problems with sample size e.g. 500 vs 501
datapre    = ft_redefinetrial(cfg, data_preprocessed{1});
cfg.toilim = [0 0.35];
datapst    = ft_redefinetrial(cfg, data_preprocessed{1});

cfg                  = [];
cfg.channel          = chantype;
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avgpre               = ft_timelockanalysis(cfg,datapre);
avgpst               = ft_timelockanalysis(cfg,datapst);

cfg                       = [];
cfg.grad                  = avg{1}.grad;  % gradiometer distances
cfg.channel               = chantype;                       
cfg.sourcemodel           = sourcemodel;
cfg.headmodel             = headmodel;     
cfg.singleshell.batchsize = 2000;
% cfg.normalize             = 'yes'; 
lf                        = ft_prepare_leadfield(cfg);

% common filter
cfg                 = [];
cfg.channel         = chantype;
cfg.method          = 'lcmv';
cfg.sourcemodel     = lf; % leadfield
cfg.headmodel       = headmodel; % volume conduction model (headmodel)
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.kappa      = 70;
% cfg.lcmv.lambda     = '5%';
% cfg.lcmv.weightnorm = 'unitnoisegain';
% cfg.lcmv.weightnorm = 'nai';
% cfg.lcmv.projectnoise = 'yes';
source1             = ft_sourceanalysis(cfg, avg{1});

cfg                    = [];
cfg.channel            = chantype;
cfg.method             = 'lcmv';
cfg.sourcemodel        = lf;      
cfg.headmodel          = headmodel;       
cfg.sourcemodel.filter = source1.avg.filter; 
sourcepre              = ft_sourceanalysis(cfg, avgpre); 
sourcepst              = ft_sourceanalysis(cfg, avgpst); 

sourceDiff         = sourcepst;
sourceDiff.avg.pow = (sourcepst.avg.pow - sourcepre.avg.pow) ./ sourcepre.avg.pow;

% interpolate the source power onto the individual MRI
mri = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_T1w-segmented.mat']));
% Keep it in SI-units
mri = ft_convert_units(mri, 'm'); 
cfg              = [];
cfg.parameter    = 'avg.pow';
cfg.voxelcoord   = 'no';
cfg.interpmethod = 'nearest';
sourceDiff_int      = ft_sourceinterpolate(cfg, sourceDiff, mri);

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
% cfg.location = [-42 -18 67];
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, sourceDiff_int);




%% Clean-Up
rmpath(fullfile('..','..'))
rmpath(fullfile('..','..','helper_functions'))