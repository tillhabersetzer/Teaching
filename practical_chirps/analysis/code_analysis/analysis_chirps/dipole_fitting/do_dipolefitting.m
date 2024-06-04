 close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the averages and
% dipolefitting
% does baseline normalization
%
% without prewhitening! -> uses only gradiometers for advanced calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject07'};

% channels of interest for sourcemodelling
coi = 'meg';

% check results in plots
check = 0;

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [4, 30];

% apply ica
ica_status = 0;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% apply baseline correction
baseline_correction_status = 1;

% downsample data
downsample_status = 0;
fs_down           = 300;

% save data
save_data = 1;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);
% addpath for ica functions
addpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);
% addpath for preprocessing function
addpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);

%% analysis

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj

%% preprocessing
config.subject                    = subjects{s};
config.files2preproc              = files2preproc;
config.bp_freq                    = bp_freq;
config.ica_status                 = ica_status;
config.ica_files                  = ica_files;
config.baseline_correction_status = baseline_correction_status;
config.downsample_status          = downsample_status;
config.fs_down                    = fs_down;

data_preprocessed = preprocess_chirps(config);
eval(subjects{s})
    
cfg = [];
avg = ft_timelockanalysis(cfg,data_preprocessed);

%% Fit a dipole model to the MEG data
%--------------------------------------------------------------------------
% template_grid = importdata(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ...
%                             'template\sourcemodel\standard_sourcemodel3d10mm.mat']);
% template_grid = ft_convert_units(template_grid,'mm');

headmodel     = importdata(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_headmodel.mat'])); % mm
sourcemodel   = importdata(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_grid.mat'])); % mm
mri_segmented = importdata(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_mri_segmented.mat'])); % mm

% headmodel   = ft_convert_units(headmodel, data_preprocessed.grad.unit);
% sourcemodel = ft_convert_units(sourcemodel, data_preprocessed.grad.unit);

% check coordinate systems
%-------------------------
if check
    dataset = fullfile(subjectdata.rawdatadir,'fixer_olsa_tsss_mc.fif');
    shape   = ft_read_headshape(dataset,'unit','mm'); 

    figure
    ft_plot_headmodel(headmodel,'facealpha',0.5);
    hold on
    ft_plot_sens(ft_convert_units(data_preprocessed.grad,'mm'), 'style', '*b');
    ft_plot_headshape(shape);
    ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:), 'maskstyle', 'opacity', 'facecolor', 'black', ...
                 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);
    view ([0 -90 0])
    title(subjects{s})
end

% prepare leadfield - does not work if one select different channels
% afterwars
%------------------
% cfg                       = [];
% cfg.grad                  = data_preprocessed.grad;    % sensor information
% cfg.channel               = coi;                       % the used channels
% cfg.sourcemodel           = sourcemodel;               % source points
% cfg.headmodel             = headmodel;                 % volume conduction model
% cfg.singleshell.batchsize = 5000;                      % speeds up the computation
% leadfield                 = ft_prepare_leadfield(cfg); % NOTE: input of the whitened data ensures the correct sensor definition to be used.

% symmetrical grid search
%------------------------
cfg             = [];
cfg.latency     = [0.07 0.14];
cfg.numdipoles  = 2;
cfg.symmetry    = 'x';
cfg.sourcemodel = sourcemodel;
cfg.gridsearch  = 'yes';
cfg.headmodel   = headmodel; 
cfg.senstype    = 'meg';
cfg.channel     = 'megplanar';
source_planar   = ft_dipolefitting(cfg,avg);

cfg.channel     = 'megmag';
source_mag      = ft_dipolefitting(cfg,avg);

cfg.channel     = 'meg';
source_meg      = ft_dipolefitting(cfg,avg);

% check
%------
if check

figure
hold on

ft_plot_dipole(source_mag.dip.pos(1,:), mean(source_mag.dip.mom(1:3,:),2), 'color', 'r', 'unit', 'mm')
ft_plot_dipole(source_mag.dip.pos(2,:), mean(source_mag.dip.mom(4:6,:),2), 'color', 'r', 'unit', 'mm')

ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g', 'unit', 'mm')
ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g', 'unit', 'mm')

ft_plot_dipole(source_meg.dip.pos(1,:), mean(source_meg.dip.mom(1:3,:),2), 'color', 'b', 'unit', 'mm')
ft_plot_dipole(source_meg.dip.pos(2,:), mean(source_meg.dip.mom(4:6,:),2), 'color', 'b', 'unit', 'mm')

pos = mean(source_mag.dip.pos,1);
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)

%savefig([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'dipolefit_megplanar_megmag_meg_sym.fig'])
end


% use estimated positions as a new starting point
%------------------------------------------------
cfg                 = [];
cfg.latency         = [0.070 0.140];
cfg.numdipoles      = 2;
cfg.symmetry        = [];
cfg.gridsearch      = 'no';
cfg.dip.pos         = source_planar.dip.pos;
cfg.sourcemodel     = sourcemodel;
cfg.headmodel       = headmodel;
cfg.channel         = 'megplanar';
cfg.senstype        = 'meg';
source_nosym        = ft_dipolefitting(cfg,avg);

% check
%------
if check 

figure
hold on

ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g', 'unit', 'mm')
ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g', 'unit', 'mm')

ft_plot_dipole(source_nosym.dip.pos(1,:), mean(source_nosym.dip.mom(1:3,:),2), 'color', 'm',  'unit', 'mm')
ft_plot_dipole(source_nosym.dip.pos(2,:), mean(source_nosym.dip.mom(4:6,:),2), 'color', 'm',  'unit', 'mm')

pos = mean(source_planar.dip.pos,1);
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)

%savefig([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'dipolefit_megplanar_sym_nosym.fig'])
end

% estimate timecourse activity
%-----------------------------
cfg = [];
cfg.latency     = 'all';
cfg.numdipoles  = 2;
cfg.symmetry    = [];
cfg.nonlinear   = 'no';  % use a fixed position
cfg.gridsearch  = 'no';
cfg.dip.pos     = source_nosym.dip.pos;
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel;
cfg.channel     = 'megplanar';
cfg.senstype    = 'meg';
source_all      = ft_dipolefitting(cfg, avg); % estimate the amplitude and orientation

% check
%------
if check 

figure
c = {'b','g','r'};
subplot(2,1,1); title('megplanar: probably right')
hold on
%plot(source_all.time, source_all.dip.mom(1:3,:), '-','Color',{'r','y','b'})
arrayfun(@(i) plot(source_all.time,source_all.dip.mom(i,:),'-','color',c{i}),1:3)
legend({'x', 'y', 'z'});
%axis([-0.1 0.4 -40e-3 40e-3])
grid on

subplot(2,1,2); title('megplanar: probably left')
hold on
%plot(source_all.time, source_all.dip.mom(4:6,:), '-')
arrayfun(@(i) plot(source_all.time,source_all.dip.mom(i,:),'-','color',c{i-3}),4:6)
legend({'x', 'y', 'z'});
%axis([-0.1 0.4 -40e-3 40e-3])
grid on

%savefig([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'dipolefit_timecourse_megplanar_nosym.fig'])
end

% moving dipole model
%--------------------
cfg             = [];
cfg.model       = 'moving';  % default is rotating
cfg.latency     = [0.07 0.14]; % [0.050 0.20]
cfg.numdipoles  = 2;
cfg.gridsearch  = 'no';
cfg.dip.pos     = source_nosym.dip.pos;
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel;
cfg.channel     = 'megplanar';
cfg.senstype    = 'meg';
sources         = ft_dipolefitting(cfg, avg);


% check
%------
if check
% copy the time-varying position of the two dipoles into a single matrix for convenience.
for i=1:numel(sources.dip)
pos1(i,:) = sources.dip(i).pos(1,:);
pos2(i,:) = sources.dip(i).pos(2,:);
end

figure
hold on
plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'r.')
plot3(pos2(:,1), pos2(:,2), pos2(:,3), 'g.')
pos = (mean(pos1, 1) + mean(pos2, 1))/2;

ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off

%savefig([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'moving_dipolefit_megplanar.fig'])
end

close all;

results_dipolefitting.source            = source;       % fitted for specified latency
results_dipolefitting.source_nosym      = source_nosym; % fitted for specified latency
results_dipolefitting.source_timeseries = source_all;   % fitted for complete interval
results_dipolefitting.moving_source     = sources;      % moving dipole
save(fullfile(subjectdata.chirps_dipolefitting,files2preproc,[subjectdata.subjectname,'_results_dipolefitting.mat']),'results_dipolefitting')

end % subjects

disp('do_dipolefitting finished')

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files']);
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);
rmpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);