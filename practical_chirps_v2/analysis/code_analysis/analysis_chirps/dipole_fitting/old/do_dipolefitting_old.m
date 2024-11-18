 close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the averages and
% dipolefitting
% does baseline normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject04'};

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [1, 30];

% apply ica
ica_on = 1;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% downsample data
downsample_data = 1;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);
% addpath for ica functions
addpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);

%% initialize

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj
    % subject selected
subject = subjects{s};
eval(subject)

filenames = get_filenames(subjectdata,files2preproc);
N_files   = length(filenames);

% loop over selected files
%-------------------------
data_preprocessed = cell(1,N_files);

%parfor f = 1:N_files
for f = 1:N_files
    
    path_dataset = [subjectdata.rawdatadir filesep filenames{f} '.fif'];
    % filter continuous data to avoid edge artifacts
    cfg              = [];
    cfg.dataset      = path_dataset;
    cfg.channel      = 'meg'; 
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;            % ensure that sensors are expressed in SI units
    data             = ft_preprocessing(cfg);   
    
    %% reject earlier specified independet components
    if ica_on
        data = reject_independent_components(data,subjectdata,ica_files);
    end
    
    %% filter data
    cfg              = [];
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = bp_freq;
    cfg.dftfilter    = 'yes';        % enable notch filtering to eliminate power line noise
    cfg.dftfreq      = [50 100 150]; % set up the frequencies for notch filtering
    cfg.coilaccuracy = 0;
    data             = ft_preprocessing(cfg,data);   

    %% define trials
    cfg                     = [];
    cfg.dataset             = path_dataset;
    cfg.trialfun            = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = 2;
    cfg.trialdef.prestim    = 0.5;                  % in seconds
    cfg.trialdef.poststim   = 0.5;                  % in seconds
    cfg                     = ft_definetrial(cfg);
    trl                     = cfg.trl;
    
    %% detect bad trials
    cfg                        = [];
    cfg.trl                    = trl;
    cfg.dataset                = path_dataset;
    cfg.artfctdef.jump.channel = 'meg'; 
    [~, artifact_jump]         = ft_artifact_jump(cfg);
    
    cfg                        = [];
    cfg.trl                    = trl;
    cfg.dataset                = path_dataset;
    cfg.artfctdef.clip.channel = 'meg'; 
    [~, artifact_clip]         = ft_artifact_clip(cfg);
    
    cfg                         = [];
    cfg.trl                     = trl;
    cfg.dataset                 = path_dataset;
    cfg.artfctdef.jump.artifact = artifact_jump;
    cfg.artfctdef.clip.artifact = artifact_clip;
    cfg = ft_rejectartifact(cfg);
    trl_new = cfg.trl;
    
    %% epoch data
    cfg                = [];
    cfg.trl            = trl_new;          
    data               = ft_redefinetrial(cfg,data);
    
    %% baseline correction
    cfg                = [];
    cfg.baselinewindow = [-0.5 0];
    cfg.demean         = 'yes';
    data               = ft_preprocessing(cfg,data); 
    
    %% again trial rejection for the rest
    % better to do it after removing of the hugh ecg artifacts
    % magnetometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = 'megmag'; 
    cfg.artfctdef.threshold.range    = 2500*10^-15; % 1000 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold1]         = ft_artifact_threshold(cfg,data);
    
    % gradiometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = 'megplanar'; 
    cfg.artfctdef.threshold.range    = 2500*10^-15/(4*10^-2); % 800 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold2]         = ft_artifact_threshold(cfg,data);
    
    cfg                              = [];
    cfg.artfctdef.threshold.artifact = [artifact_threshold1;artifact_threshold2];
    data = ft_rejectartifact(cfg,data);
    
    %% downsample data
    if downsample_data
        cfg            = [];
        cfg.resamplefs = 300;
        cfg.detrend    = 'no';
        data           = ft_resampledata(cfg,data);
    end
  
    data_preprocessed{f} = data;
    
    clear data

end % files
    
%% append data
hdr                   = data_preprocessed{1}.hdr;
cfg                   = [];
cfg.keepsampleinfo    = 'no';
data_preprocessed     = ft_appenddata(cfg,data_preprocessed{:});
data_preprocessed.hdr = hdr;

%% average data
cfg      = [];
data_avg = ft_timelockanalysis(cfg,data_preprocessed);

%% combine planar gradient
cfg            = [];
cfg.method     = 'sum';
cfg.updatesens = 'yes'; % need old information for dipole fitting 'no'
data_avg_cmb   = ft_combineplanar(cfg,data_avg);

%% plot data
    
new_dir = [subjectdata.chirps_dipolefitting filesep files2preproc];
if ~exist(new_dir, 'dir')
    mkdir(new_dir)
end

% magnetometer
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,data_avg_cmb);
suptitle('avergage magnetometer')
savefig([new_dir filesep 'avg_mag.fig'])

figure
cfg          = [];
cfg.xlim     = [0.07, 0.140];
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306mag.lay';
ft_topoplotER(cfg,data_avg_cmb);
suptitle('topoplot magnetometer')
savefig([new_dir filesep 'avg_topo_mag.fig'])

figure
cfg        = [];
cfg.xlim   = [-0.1 : 0.1 : 0.5];  % Define time intervals
cfg.layout = 'neuromag306mag.lay';
ft_topoplotER(cfg,data_avg_cmb);
suptitle('topoplot series magnetometer')
savefig([new_dir filesep 'series_topo_mag.fig'])

% planar gradient
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg, data_avg_cmb);
suptitle('average combined gradiometers')
savefig([new_dir filesep 'avg_cmbplanar.fig'])

figure
cfg          = [];
cfg.xlim     = [0.07, 0.140];
cfg.colorbar = 'yes';
cfg.layout   = 'neuromag306cmb.lay';
ft_topoplotER(cfg,data_avg_cmb);
suptitle('topoplot combined gradiometers')
savefig([new_dir filesep 'avg_topo_cmbplanar.fig'])

figure
cfg        = [];
cfg.xlim   = [-0.1 : 0.1 : 0.5];  % Define time intervals
cfg.layout = 'neuromag306cmb.lay';
ft_topoplotER(cfg,data_avg_cmb);
suptitle('topoplot series combined gradiometers')
savefig([new_dir filesep 'series_topo_cmbplanar.fig'])

close all

%% Fit a dipole model to the MEG data
%--------------------------------------------------------------------------
template_grid = importdata(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ...
                            'template\sourcemodel\standard_sourcemodel3d10mm.mat']);
template_grid = ft_convert_units(template_grid,'mm');

% load atlas and create binary mask (restrict your template grid further with a mask)
atlas = ft_read_atlas('Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\template\atlas\aal\ROI_MNI_V4.nii');

cfg            = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;
cfg.inputcoord = 'mni';
mask           = ft_volumelookup(cfg,template_grid);

% create temporary mask according to the atlas entries
tmp         = repmat(template_grid.inside,1,1);
tmp(tmp==1) = 0;
tmp(mask)   = 1;

% define inside locations according to the atlas based mask
template_grid.inside = tmp;

% plot the atlas based grid
figure; 
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

% Inverse-warp the subject specific grid to the atlas based template grid (make individual subject's grid)
mri_segmented = importdata([subjectdata.headmodel filesep 'mri_segmented.mat']);

cfg           = [];
cfg.warpmni   = 'yes';
cfg.template  = template_grid;
cfg.nonlinear = 'yes';
cfg.mri       = mri_segmented;
cfg.unit      = 'mm';
sourcemodel   = ft_prepare_sourcemodel(cfg);

headmodel_meg    = importdata([subjectdata.headmodel filesep 'headmodel_meg.mat']);

% symmetrical grid search
cfg             = [];
cfg.latency     = [0.07 0.14];
cfg.numdipoles  = 2;
cfg.symmetry    = 'x';
cfg.sourcemodel = sourcemodel;
cfg.gridsearch  = 'yes';
cfg.headmodel   = headmodel_meg;
cfg.senstype    = 'meg';
cfg.channel     = 'megplanar';
source_planar   = ft_dipolefitting(cfg, data_avg);

cfg.channel    = 'megmag';
source_mag     = ft_dipolefitting(cfg, data_avg);

cfg.channel    = 'meg';
source_meg     = ft_dipolefitting(cfg, data_avg);

% check
%------
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

savefig([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'dipolefit_megplanar_megmag_meg_sym.fig'])

% use estimated positions as a new starting point
%------------------------------------------------
% - only for magnetometers
cfg                 = [];
cfg.latency         = [0.070 0.140];
cfg.numdipoles      = 2;
cfg.symmetry        = [];
cfg.gridsearch      = 'no';
cfg.dip.pos         = source_planar.dip.pos;
cfg.sourcemodel     = sourcemodel;
cfg.headmodel       = headmodel_meg;
cfg.channel         = 'megplanar';
cfg.senstype        = 'meg';
source_planar_nosym = ft_dipolefitting(cfg, data_avg);

% check
%------
figure
hold on

ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g', 'unit', 'mm')
ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g', 'unit', 'mm')

ft_plot_dipole(source_planar_nosym.dip.pos(1,:), mean(source_planar_nosym.dip.mom(1:3,:),2), 'color', 'm',  'unit', 'mm')
ft_plot_dipole(source_planar_nosym.dip.pos(2,:), mean(source_planar_nosym.dip.mom(4:6,:),2), 'color', 'm',  'unit', 'mm')

pos = mean(source_planar.dip.pos,1);
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)

savefig([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'dipolefit_megplanar_sym_nosym.fig'])

% estimate timecourse activity
%-----------------------------
cfg = [];
cfg.latency     = 'all';
cfg.numdipoles  = 2;
cfg.symmetry    = [];
cfg.nonlinear   = 'no';  % use a fixed position
cfg.gridsearch  = 'no';
cfg.dip.pos     = source_planar_nosym.dip.pos;
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel_meg;
cfg.channel     = 'megplanar';
cfg.senstype    = 'meg';
source_all      = ft_dipolefitting(cfg, data_avg); % estimate the amplitude and orientation

% check
%------
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

savefig([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'dipolefit_timecourse_megplanar_nosym.fig'])

% moving dipole model
%--------------------
cfg             = [];
cfg.model       = 'moving';  % default is rotating
cfg.latency     = [0.07 0.14]; % [0.050 0.20]
cfg.numdipoles  = 2;
cfg.gridsearch  = 'no';
cfg.dip.pos     = source_planar_nosym.dip.pos;
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel_meg;
cfg.channel     = 'megplanar';
cfg.senstype    = 'meg';
source          = ft_dipolefitting(cfg, data_avg);

% copy the time-varying position of the two dipoles into a single matrix for convenience.
for i=1:numel(source.dip)
pos1(i,:) = source.dip(i).pos(1,:);
pos2(i,:) = source.dip(i).pos(2,:);
end

% check
%------
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

savefig([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'moving_dipolefit_megplanar.fig'])

close all;

results_dipolefitting.source_planar       = source_planar;       % fitted for specified latency
results_dipolefitting.source_planar_nosym = source_planar_nosym; % fitted for specified latency
results_dipolefitting.source_timeseries   = source_all;          % fitted for complete interval
results_dipolefitting.moving_source       = source;              % moving dipole
save([subjectdata.chirps_dipolefitting filesep files2preproc filesep 'results_dipolefitting.mat'],'results_dipolefitting')

end % subjects

disp('do_dipolefitting finished')

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])

 
% have a look
% before ica
% cfg                              = [];
% cfg.channel                      = 'megplanar';
% cfg.dataset                      = [subjectdata.rawdatadir filesep filenames{1} '.fif'];
% cfg.blocksize                    = 100; % pulses are in the first 100sec
% cfg.artfctdef.threshold.artifact = artifact_threshold1; % use thresholds for correct file
% cfg.plotevents                   = 'no';
% ft_databrowser(cfg)

% after ica
% cfg                         = [];
% cfg.channel                 = 'megplanar';
% cfg.continuous              = 'yes';
% cfg.blocksize               = 100;
% cfg.artfctdef.threshold.artifact = artifact_threshold1;
% cfg.plotevents              = 'no';
% ft_databrowser(cfg,data_preprocessed) % load preprocessed data with applied ica

