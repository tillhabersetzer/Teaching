%--------------------------------------------------------------------------
% Check out: https://www.fieldtriptoolbox.org/workshop/natmeg/dipolefitting/
%
% This script performas dipolfitting with a loop over all subjects
%
%--------------------------------------------------------------------------

close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','..','subjectdata'))
eval('main_settings')

addpath(fullfile('..','..','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
subjects = 4;

% Choose which type of evoked field (erf) to fit
erf_type = 'N19mP30m'; % type for practical
% erf_type = 'N100m'; 

% choose sensortype
% chantype = {'meg','eeg'};
chantype = {'eeg'};
% chantype = {'meg'};

% type of sourcemodel
% sourcemodel_type = 'vol';
sourcemodel_type = 'surf';

% Switch between erf-types
%-------------------------
switch erf_type
    case 'N19mP30m'
        timewin = [0.02,0.05]; % 20-50ms covering N19m-P30m-Peak2Peak interval
    case 'N100m'
        timewin = [0.095,0.125]; % N100m peak
end 

Nsub  = length(subjects);
Nsens = length(chantype);

%% Preprocess data
%--------------------------------------------------------------------------

for subidx=subjects % loop over subjects

    subject  = ['sub-',num2str(subidx,'%02d')];

    % Load sourcemodel
    %----------------- 
    % Proceed in SI-units. Leads apparantly to correct units for dipole moments
    % in Am
    
    switch sourcemodel_type
        case 'vol'
            sourcemodel = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','sourcemodel',[subject,'_sourcemodel-volumetric.mat']));
        case 'surf'
            sourcemodel = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','sourcemodel',[subject,'_sourcemodel-corticalsheet4k.mat']));
    end

    sourcemodel = ft_convert_units(sourcemodel, 'm'); 

    source_vec      = cell(1,3); % 3 conditions
    source_sca_mean = cell(1,3); 
    source_sca_svd  = cell(1,3);  

    for sidx=1:Nsens % loop over channel types

        % Load data
        %----------
        data       = importdata(fullfile(settings.path2bids,'derivatives',subject,'sensorlevel',[subject,'_erf-',erf_type,'_',chantype{sidx},'.mat']));
        avg        = data.avg;
        avg_pooled = data.avg_pooled;
        clear data
    
        %% Fit initial dipole model to pooled conditions (clicks+downchirps)
        %----------------------------------------------------------------------
        
        switch chantype{sidx}
            case {'meg','megmag','megplanar'}
                headmodel = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_headmodel-singleshell.mat']));
                headmodel = ft_convert_units(headmodel, 'm');  
                senstype  = 'meg';
            case 'eeg'
                headmodel = importdata(fullfile(settings.path2bids,'derivatives',subject,'forward_modelling','headmodel',[subject,'_headmodel-bem.mat']));
                headmodel = ft_convert_units(headmodel, 'm');  
                senstype  = 'eeg';

                % need to install Openmeeg in case method 'openmeeg' is used
                setenv('PATH', fullfile(settings.path2openmeeg,'bin'))
                setenv('LD_LIBRARY_PATH', fullfile(settings.path2openmeeg,'lib'))

                % Construction of symmetric leadfield
                %------------------------------------
                % use precomputed leadfields for eeg to speed up the
                % computation
                cfg             = [];
                cfg.headmodel   = headmodel;
                cfg.sourcemodel = sourcemodel;
                cfg.channel     = 'eeg'; 
                cfg.elec        = avg_pooled.elec;
                leadfield_sym   = compute_sym_leadfield(cfg);
        end

        % Apply Sphering if magnetometers and gradiometers are used
        % together
        %----------------------------------------------------------
        sphering = strcmp(chantype{sidx},'meg');
        if sphering
           filenames   = {'noisepre','noisepost'};
           % check if file exists
           noisepath1 = fullfile(settings.path2bids,'derivatives',subject,'tsss',[subject,'_task-',filenames{1},'_meg-raw_tsss.fif']);
           noisepath2 = fullfile(settings.path2bids,'derivatives',subject,'tsss',[subject,'_task-',filenames{2},'_meg-raw_tsss.fif']);
           datapath   = [];
           if isfile(noisepath1)
               datapath = {noisepath2};
           end
           if isfile(noisepath2)
               datapath{end+1} = noisepath2;
           end
           avg_noise = compute_noisecov(datapath,'meg');
        end

        % Symmetrical gridsearch
        %-----------------------
        cfg             = [];
        cfg.latency     = timewin;
        cfg.numdipoles  = 2; % fit 2 dipoles
        cfg.symmetry    = 'x'; % expect activity in both auditory cortices, fit symmetric dipoles
        if strcmp(chantype{sidx},'eeg')
            cfg.sourcemodel = leadfield_sym; % use procomputed leadfields for eeg
        else
            cfg.sourcemodel = sourcemodel;
        end
        cfg.gridsearch  = 'yes';
        cfg.headmodel   = headmodel;
        cfg.senstype    = senstype;
        cfg.channel     = chantype{sidx}; 
        
        % Optionally, you can include a noise covariance structure to sphere the data (is useful when using both magnetometers and gradiometers to fit your dipole)
        if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
        source_pooled = ft_dipolefitting(cfg, avg_pooled);

        % Now that we have a better starting point for the dipole fit, we can release the symmetry contstraint
        %-----------------------------------------------------------------------------------------------------
        cfg            = [];
        cfg.latency    = timewin;
        cfg.numdipoles = 2;
        cfg.symmetry   = []; % no symmetry constraint
        cfg.gridsearch = 'no';
        cfg.nonlinear  = 'yes'; % non-linear search for optimal dipole parameters
        cfg.dip.pos    = source_pooled.dip.pos;
        cfg.sourcemodel.unit = sourcemodel.unit; % defines units of dipole positions
        cfg.headmodel  = headmodel;
        cfg.channel    = chantype{sidx}; 
        cfg.senstype   = senstype;
        % Optionally, you can include a noise covariance structure to sphere the data (is useful when using both magnetometers and gradiometers to fit your dipole)
        if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
        source_pooled_nosym = ft_dipolefitting(cfg, avg_pooled);

        %% Reconstruct time courses of all conditions with different approaches
        %----------------------------------------------------------------------
        
        % 1. Approach: Fit freely oriented dipole (fixed positions, loose orientation)
        %-----------------------------------------------------------------------------
        cfg              = [];
        cfg.latency      = 'all';
        cfg.numdipoles   = 2;
        cfg.symmetry     = [];
        cfg.nonlinear    = 'no';  % use a fixed position
        cfg.gridsearch   = 'no';
        cfg.dip.pos      = source_pooled_nosym.dip.pos; % use estimated dipole positions from pooled data as reference
        cfg.sourcemodel.unit = sourcemodel.unit; % defines units of dipole positions
        cfg.headmodel    = headmodel;
        cfg.channel      = chantype{sidx}; 
        cfg.senstype     = senstype;
        if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
        source_vec{1} = ft_dipolefitting(cfg, avg{1}); % estimate the amplitude and orientation
        source_vec{2} = ft_dipolefitting(cfg, avg{2}); % estimate the amplitude and orientation
        source_vec{3} = ft_dipolefitting(cfg, avg{3}); % estimate the amplitude and orientation

    
        % 2. Approach: Fit fixed oriented dipole (fixed positions, fixed orientation)
        %----------------------------------------------------------------------------
        % dipole orientations are estimated based on previously computed
        % dipolmoments (xyz-components)
        % mean orientations is used -> returns single dipolar timecourse
        cfg           = [];
        cfg.headmodel = headmodel;
        cfg.channel   = chantype{sidx};     
        cfg.dippos    = source_pooled_nosym.dip.pos;
        cfg.dipmom    = source_pooled_nosym.dip.mom;
        cfg.unit      = sourcemodel.unit;
        if sphering; cfg.noisecov = avg_noise.cov; end  
        
        source_sca_mean{1} = constrained_dipolfitting(cfg,avg{1});
        source_sca_mean{2} = constrained_dipolfitting(cfg,avg{2});
        source_sca_mean{3} = constrained_dipolfitting(cfg,avg{3});
    
        % 3. Approach: Fit fixed orienteddipole (fixed positions, fixed orientation)
        %---------------------------------------------------------------------------
        % dipole orientations are estimated based on previously computed
        % dipolmoments (xyz-components) 
        % maximum variance orientation via SVD is used to return single dipolar timecourse
        
        % Calculate dipole moment orientations
        %-------------------------------------
        moments   = source_pooled_nosym.dip.mom(1:3,:);
        [u, ~, ~] = svd(moments, 'econ'); 
        ori1      = u(:,1); % already normalized
        moments   = source_pooled_nosym.dip.mom(4:6,:);
        [u, ~, ~] = svd(moments, 'econ'); 
        ori2      = u(:,1); % already normalized
        
        cfg           = [];
        cfg.dippos    = source_pooled_nosym.dip.pos;
        cfg.headmodel = headmodel;
        cfg.channel   = chantype{sidx};
        cfg.ori       = [ori1,ori2];
        cfg.unit      = sourcemodel.unit;   

        source_sca_svd{1} = constrained_dipolfitting(cfg,avg{1});
        source_sca_svd{2} = constrained_dipolfitting(cfg,avg{2});
        source_sca_svd{3} = constrained_dipolfitting(cfg,avg{3});


        %% Create structure for saving data
        %----------------------------------

        % Save reconstructed waveforms and dipole positions
        %--------------------------------------------------
        dir2save = fullfile(settings.path2bids,'derivatives',subject,'sourcelevel');
        if ~exist(dir2save, 'dir')
            mkdir(dir2save)
        end
    
        data                      = [];
        % dipole positions
        data.dippos.pos_sym       = source_pooled.dip.pos;
        data.dippos.pos_nosym     = source_pooled_nosym.dip.pos;
        % dipole moments
        data.source_pooled        = source_pooled;
        data.source_pooled_nosym  = source_pooled_nosym;
        % vector dipol moments of conditions
        data.source_vec           = source_vec;
        % scalar dipol moments of conditions
        data.source_sca_mean      = source_sca_mean;
        % scalar dipol moments of conditions
        data.source_sca_svd       = source_sca_svd;
        % metadata
        data.conditions           = {'clicks','upchirps','downchirps'};
        data.sourcemodel_type     = sourcemodel_type;
        
        save(fullfile(dir2save,[subject,'_dipolefits-',erf_type,'_',chantype{sidx},'.mat']),'-struct','data');
    
        clear source_pooled source_pooled_nosym source_vec source_sca_mean source_sca_svd data

    end % loop over sensor types

end % end loop over subjects
