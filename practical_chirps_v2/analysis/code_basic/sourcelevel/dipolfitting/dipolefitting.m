%--------------------------------------------------------------------------
% Check out: https://www.fieldtriptoolbox.org/workshop/natmeg/dipolefitting/
%
% This script performs dipolfitting with a loop over all subjects
%--------------------------------------------------------------------------

close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','..'))
eval('main_settings')

addpath(fullfile('..','..','helper_functions'))

% check whether maxfiltered data should be analyzed
maxfilter = settings.maxfilter;

%% Script settings
%--------------------------------------------------------------------------
subjects = [1,2];

% Switch between erf-types
%-------------------------
timewin = [0.02 0.05]; % 20-50ms covering N19m-P30m-Peak2Peak interval

chantype = 'meg'; % 'megmag','megplanar','meg'

%% Preprocess data
%--------------------------------------------------------------------------

for subidx=subjects % loop over subjects

    subject  = ['sub-',num2str(subidx,'%02d')];

    % Load sourcemodel and headmodel
    %-------------------------------
    % Proceed in SI-units. Leads apparantly to correct units for dipole moments in Am

    sourcemodel = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_sourcemodel-volumetric.mat']));
    sourcemodel = ft_convert_units(sourcemodel, 'm'); 

    headmodel = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_headmodel-singleshell.mat']));
    headmodel = ft_convert_units(headmodel, 'm');  

    source_vec          = cell(1,3); % 3 conditions
    source_sca_mean     = cell(1,3); 
    source_sca_svd      = cell(1,3);  

    % Load data
    %----------
    data       = importdata(fullfile(settings.path2project,'derivatives',subject,'sensorlevel',[subject,'_erf-N19mP30m.mat']));
    avg        = data.avg;
    avg_pooled = data.avg_pooled;
    clear data
    
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

%           imagesc(avg_noise.cov)
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

    %% Reconstruct time courses of all conditions with different approaches
    %----------------------------------------------------------------------
    
    % 1. Approach: Fit freely oriented dipole (fixed positions, loose orientation)
    %-----------------------------------------------------------------------------
    cfg                  = [];
    cfg.latency          = 'all';
    cfg.numdipoles       = 2;
    cfg.symmetry         = [];
    cfg.nonlinear        = 'no';  % use a fixed position
    cfg.gridsearch       = 'no';
    cfg.dip.pos          = source_pooled_nosym.dip.pos; % use estimated dipole positions from pooled data as reference
    cfg.sourcemodel.unit = sourcemodel.unit; % defines units of dipole positions
    cfg.headmodel        = headmodel;
    cfg.channel          = chantype; 
    cfg.senstype         = 'meg';
    
    if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
    source_vec{1} = ft_dipolefitting(cfg, avg{1}); % estimation of amplitude and orientation
    source_vec{2} = ft_dipolefitting(cfg, avg{2}); 
    source_vec{3} = ft_dipolefitting(cfg, avg{3}); 

    % 2. Approach: Fit fixed oriented dipole (fixed positions, fixed orientation)
    %----------------------------------------------------------------------------
    % dipole orientations are estimated based on previously computed
    % dipolmoments (xyz-components)
    % mean orientation is used to return single dipolar timecourse
    cfg           = [];
    cfg.headmodel = headmodel;
    cfg.channel   = chantype;     
    cfg.dippos    = source_pooled_nosym.dip.pos;
    cfg.dipmom    = source_pooled_nosym.dip.mom;
    cfg.diptime   = source_pooled_nosym.time;
    cfg.unit      = sourcemodel.unit;
    cfg.ori       = 'mean';
    cfg.timewin   = timewin;
    
    if sphering; cfg.noisecov = avg_noise.cov; end 
    source_sca_mean{1} = constrained_dipolfitting(cfg,avg{1});
    source_sca_mean{2} = constrained_dipolfitting(cfg,avg{2});
    source_sca_mean{3} = constrained_dipolfitting(cfg,avg{3});

    % 3. Approach: Fit fixed oriented dipole (fixed positions, fixed orientation)
    %----------------------------------------------------------------------------
    % dipole orientations are estimated based on previously computed
    % dipolmoments (xyz-components) 
    % maximum variance orientation via SVD is used to return single dipolar timecourse
    cfg           = [];
    cfg.headmodel = headmodel;
    cfg.channel   = chantype;     
    cfg.dippos    = source_pooled_nosym.dip.pos;
    cfg.dipmom    = source_pooled_nosym.dip.mom;
    cfg.diptime   = source_pooled_nosym.time;
    cfg.unit      = sourcemodel.unit;
    cfg.ori       = 'svd';
    cfg.timewin   = timewin;
    
    if sphering; cfg.noisecov = avg_noise.cov; end 
    source_sca_svd{1} = constrained_dipolfitting(cfg,avg{1});
    source_sca_svd{2} = constrained_dipolfitting(cfg,avg{2});
    source_sca_svd{3} = constrained_dipolfitting(cfg,avg{3});

    %% Create structure for saving data
    %----------------------------------

    % Save reconstructed waveforms and dipole positions
    %--------------------------------------------------
    dir2save = fullfile(settings.path2project,'derivatives',subject,'sourcelevel');
    if ~exist(dir2save, 'dir')
        mkdir(dir2save)
    end

    data                     = [];
    data.chantype            = chantype;
    data.conditions          = {'clicks','upchirps','downchirps'};
    % dipole moments
    data.source_pooled_sym   = source_pooled_sym;
    data.source_pooled_nosym = source_pooled_nosym;
    % vector dipol moments of conditions
    data.source_vec          = source_vec;
    % scalar dipol moments of conditions
    data.source_sca_mean     = source_sca_mean;
    % scalar dipol moments of conditions
    data.source_sca_svd      = source_sca_svd;
    
    save(fullfile(dir2save,[subject,'_dipolefits-N19mP30m.mat']),'-struct','data');

    clear source_pooled_sym source_pooled_nosym source_vec source_sca_mean source_sca_svd data

end % loop over subjects

%% Clean-Up
rmpath(fullfile('..','..'))
rmpath(fullfile('..','..','helper_functions'))