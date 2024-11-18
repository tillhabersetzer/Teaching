close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the avarages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject02'};

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [0.5, 150];

% apply ica
ica_on = 1;

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

filenames       = get_filenames(subjectdata,files2preproc);
F               = length(filenames);
filename4saving = [filenames,{files2preproc}];
%     path2save = subjectdata.ica_stories;
N_files         = F;

% loop over selected files
%-------------------------
data_preprocessed = cell(1,N_files);

for f = 1:N_files
    path_dataset = [subjectdata.rawdatadir filesep filenames{f} '.fif'];
    %% preprocess data
    cfg              = [];
    cfg.dataset      = path_dataset;
    % channels of interest and exclude bad channels
    cfg.channel      = [{'meg'},strcat('-',subjectdata.badchannels)]; 
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;
    data_cont        = ft_preprocessing(cfg);
    
    %% reject earlier specified independet components
    if ica_on
        data_cont = reject_independent_components(data_cont,subjectdata,files2preproc);
    end
    
    %% filter data    
    cfg            = [];
    cfg.channel    = 'meg';        
    cfg.detrend    = 'yes';
    cfg.demean     = 'yes';         
    cfg.continuous = 'yes';
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = bp_freq;      
    cfg.dftfilter  = 'yes';        
    cfg.dftfreq    = [50 100 150];  
    data_filt      = ft_preprocessing(cfg, data_cont);

    %% define trials
    cfg                     = [];
    cfg.dataset             = path_dataset;
    cfg.trialfun            = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype  = 'STI101';
    cfg.trialdef.eventvalue = 2;
    cfg.trialdef.prestim    = 0.05;                  % in seconds
    cfg.trialdef.poststim   = 0.3;                   % in seconds
    cfg                     = ft_definetrial(cfg);
    trl                     = cfg.trl;
    
    %% detect bad trials
    cfg                        = [];
    cfg.trl                    = trl;
    cfg.dataset                = path_dataset;
    cfg.artfctdef.jump.channel = [{'meg'},strcat('-',subjectdata.badchannels)]; 
    [~, artifact_jump]         = ft_artifact_jump(cfg);
    
    cfg                        = [];
    cfg.trl                    = trl;
    cfg.dataset                = path_dataset;
    cfg.artfctdef.clip.channel = [{'meg'},strcat('-',subjectdata.badchannels)]; 
    [~, artifact_clip]         = ft_artifact_clip(cfg);
    
    cfg                         = [];
    cfg.trl                     = trl;
    cfg.dataset                 = path_dataset;
    cfg.artfctdef.jump.artifact = artifact_jump;
    cfg.artfctdef.clip.artifact = artifact_clip;
    cfg = ft_rejectartifact(cfg);
    trl_new = cfg.trl;
     
    %% apply trial definition
    cfg     = [];
    cfg.trl = trl_new;
    data    = ft_redefinetrial(cfg, data_filt);
    
    %% again trial rejection for the rest
    % better to do it after removing of the hugh ecg artifacts
    % magnetometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = {'megmag'}; 
    cfg.artfctdef.threshold.range    = 2000*10^-15; % 1000 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold1]         = ft_artifact_threshold(cfg,data);
    
    % gradiometer
    cfg                              = [];
    cfg.artfctdef.threshold.channel  = {'megplanar'}; 
    cfg.artfctdef.threshold.range    = 2000*10^-15/(4*10^-2); % 800 fT (Stefan,Rupp)
    cfg.artfctdef.threshold.bpfilter = 'no';
    [~, artifact_threshold2]         = ft_artifact_threshold(cfg,data);
    
    cfg                              = [];
    cfg.artfctdef.threshold.artifact = [artifact_threshold1;artifact_threshold2];
    data = ft_rejectartifact(cfg,data);
    
    %% downsample data
    if downsample_data
        cfg            = [];
        cfg.resamplefs = 400;
        cfg.detrend    = 'no';
        data           = ft_resampledata(cfg,data);
    end
  
    data_preprocessed{f} = data;

end % files

%% combine data
cfg                = [];
cfg.keepsampleinfo = 'no';
data_all           = ft_appenddata(cfg,data_preprocessed{:});

% average sensor array
sens          = cellfun(@(v) v.grad,data_preprocessed);
[asens, afid] = ft_average_sens(sens);
data_all.grad = asens;

% add combined file
data_preprocessed{N_files+1} = data_all;
N_files                      = length(data_preprocessed);

%% average data
data_avg = cell(size(data_preprocessed));
cfg      = [];

for n = 1:N_files
    data_avg{n} = ft_timelockanalysis(cfg,data_preprocessed{n});
end
%% combine planar gradient
data_avg_cmb   = cell(size(data_preprocessed));
cfg            = [];
cfg.method     = 'sum';
cfg.updatesens = 'no'; % need old information for dipole fitting

for n = 1:N_files
    data_avg_cmb{n} = ft_combineplanar(cfg,data_avg{n});
end

%% plot data
for n = 1:N_files
    
    new_dir = [subjectdata.chirps_dipolefitting filesep filename4saving{n}];
    if ~exist(new_dir, 'dir')
        mkdir(new_dir)
    end
    
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306mag.lay';
    cfg.xlim       = [0, 0.2];
    ft_multiplotER(cfg,data_avg_cmb{n});
    suptitle('avergage magnetometer')
    savefig([new_dir filesep 'avg_mag.fig'])

    figure
    cfg          = [];
    cfg.xlim     = [0.1 0.3];
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306mag.lay';
    ft_topoplotER(cfg,data_avg_cmb{n});
    suptitle('topoplot magnetometer')
    savefig([new_dir filesep 'avg_topo_mag.fig'])
    
    figure
    cfg        = [];
    cfg.xlim   = [0 : 0.1 : 0.4];  % Define 12 time intervals
    cfg.layout = 'neuromag306mag.lay';
    ft_topoplotER(cfg,data_avg_cmb{n});
    suptitle('topoplot magnetometer')
    savefig([new_dir filesep 'series_topo_mag.fig'])

    % planar gradient
    figure
    cfg            = [];
    cfg.showlabels = 'yes';
    cfg.fontsize   = 6;
    cfg.layout     = 'neuromag306cmb.lay';
    ft_multiplotER(cfg, data_avg_cmb{n});
    suptitle('average combined gradiometers')
    savefig([new_dir filesep 'avg_cmbplanar.fig'])
    
    figure
    cfg          = [];
    cfg.xlim     = [0.1 0.3];
    cfg.colorbar = 'yes';
    cfg.layout   = 'neuromag306cmb.lay';
    ft_topoplotER(cfg,data_avg_cmb{n});
    suptitle('topoplot magnetometer')
    savefig([new_dir filesep 'avg_topo_cmbplanar.fig'])
    
    figure
    cfg        = [];
    cfg.xlim   = [0 : 0.1 : 0.4];  % Define 12 time intervals
    cfg.layout = 'neuromag306cmb.lay';
    ft_topoplotER(cfg,data_avg_cmb{n});
    suptitle('topoplot combined gradiometers')
    savefig([new_dir filesep 'series_topo_cmbplanar.fig'])
    
    close all
end

%% Fit a dipole model to the MEG data

% senstype = ft_senstype(data_avg.label)';
% [channel] = ft_channelselection('megplanar',data_avg.label,senstype);

% change senstype to combined planar gradiometers and test dipole fitting
% again

headmodel_meg    = importdata([subjectdata.headmodel filesep 'headmodel_meg.mat']);
mri_segmented    = importdata([subjectdata.headmodel filesep 'mri_segmented.mat']);
mri_segmented_cm = ft_convert_units(mri_segmented, 'cm');

for n = 1:N_files

% symmetrical grid search
cfg            = [];
cfg.latency    = [0.060 0.110];
cfg.numdipoles = 2;
cfg.symmetry   = 'x';
cfg.resolution = 1;
cfg.sourcemodel.unit = 'cm';
%cfg.unit       = 'cm';
cfg.gridsearch = 'yes';
cfg.headmodel  = headmodel_meg;
cfg.senstype   = 'meg';
cfg.channel    = 'megplanar';
source_planar  = ft_dipolefitting(cfg, data_avg{n});

cfg.channel    = 'megmag';
source_mag     = ft_dipolefitting(cfg, data_avg{n});

cfg.channel    = 'meg';
source_meg     = ft_dipolefitting(cfg, data_avg{n});


% check
%------
% figure
% hold on
% 
% ft_plot_dipole(source_mag.dip.pos(1,:), mean(source_mag.dip.mom(1:3,:),2), 'color', 'r')
% ft_plot_dipole(source_mag.dip.pos(2,:), mean(source_mag.dip.mom(4:6,:),2), 'color', 'r')
% 
% ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g')
% ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g')
% 
% ft_plot_dipole(source_meg.dip.pos(1,:), mean(source_meg.dip.mom(1:3,:),2), 'color', 'b')
% ft_plot_dipole(source_meg.dip.pos(2,:), mean(source_meg.dip.mom(4:6,:),2), 'color', 'b')
% 
% 
% pos = mean(source_mag.dip.pos,1);
% ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
% ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
% ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)
% 
% ft_plot_crosshair(pos, 'color', [1 1 1]/2);
% axis tight
% axis off
% view(12, -10)

% use estimated positions as a new starting point
%------------------------------------------------
cfg            = [];
cfg.latency    = [0.060 0.110];
cfg.numdipoles = 2;
cfg.symmetry   = [];
cfg.gridsearch = 'no';
cfg.dip.pos    = source_planar.dip.pos;
cfg.sourcemodel.unit = 'cm';
cfg.headmodel  = headmodel_meg;
cfg.channel    = 'megplanar';
cfg.senstype   = 'meg';
source_planar_nosym = ft_dipolefitting(cfg, data_avg{n});

% check
%------
figure
hold on

ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g')
ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g')

ft_plot_dipole(source_planar_nosym.dip.pos(1,:), mean(source_planar_nosym.dip.mom(1:3,:),2), 'color', 'm')
ft_plot_dipole(source_planar_nosym.dip.pos(2,:), mean(source_planar_nosym.dip.mom(4:6,:),2), 'color', 'm')

pos = mean(source_planar.dip.pos,1);
ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)

savefig([subjectdata.chirps_dipolefitting filesep filename4saving{n} filesep 'dipolefit_megplanar.fig'])

% estimate timecourse activity
%-----------------------------
cfg = [];
cfg.latency    = 'all';
cfg.numdipoles = 2;
cfg.symmetry   = [];
cfg.nonlinear  = 'no';  % use a fixed position
cfg.gridsearch = 'no';
cfg.dip.pos    = source_planar_nosym.dip.pos;
cfg.sourcemodel.unit = 'cm';
cfg.headmodel  = headmodel_meg;
cfg.channel    = 'megplanar';
cfg.senstype   = 'meg';
source_all     = ft_dipolefitting(cfg, data_avg{n}); % estimate the amplitude and orientation

% check
%------
figure
subplot(2,1,1); title('megplanar: left')
hold on
plot(source_all.time, source_all.dip.mom(1:3,:), '-')
legend({'x', 'y', 'z'});
%axis([-0.1 0.4 -40e-3 40e-3])
grid on

subplot(2,1,2); title('megplanar: right')
hold on
plot(source_all.time, source_all.dip.mom(4:6,:), '-')
legend({'x', 'y', 'z'});
%axis([-0.1 0.4 -40e-3 40e-3])
grid on

savefig([subjectdata.chirps_dipolefitting filesep filename4saving{n} filesep 'dipolefit_timecourse_megplanar.fig'])

% moving dipole model
%--------------------
cfg            = [];
cfg.model      = 'moving'; % default is rotating
cfg.latency    = [0.050 0.20];
cfg.numdipoles = 2;
cfg.gridsearch = 'no';
cfg.dip.pos    = source_planar_nosym.dip.pos;
cfg.sourcemodel.unit = 'cm';
cfg.headmodel  = headmodel_meg;
cfg.channel    = 'megplanar';
cfg.senstype   = 'meg';
source         = ft_dipolefitting(cfg, data_avg{n});

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

ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off

savefig([subjectdata.chirps_dipolefitting filesep filename4saving{n} filesep 'moving_dipolefit_megplanar.fig'])

close all;

results_dipolefitting.source_planar_nosym = source_planar_nosym; % fitted for specified latency
results_dipolefitting.source_timeseries   = source_all;          % fitted for complete interval
results_dipolefitting.moving_source       = source;              % moving dipole
save([subjectdata.chirps_dipolefitting filesep filename4saving{n} filesep 'results_dipolefitting.mat'],'results_dipolefitting')

end

end % subjects

disp('do_dipolefitting finished')

% % have a look
% cfg                         = [];
% cfg.channel                 = {'megmag'}; % components to be plotted
% %cfg.blocksize               = 100;
% cfg.plotevents              = 'no';
% ft_databrowser(cfg,data_preprocessed{4});
% 
% have a look
% before ica
% cfg                         = [];
% cfg.channel                 = [{'megplanar'},strcat('-',subjectdata.badchannels)];
% cfg.dataset                 = [subjectdata.rawdatadir filesep filenames{1} '.fif'];
% cfg.blocksize               = 100;
% cfg.artfctdef.threshold.artifact = artifact_threshold1;
% cfg.plotevents              = 'no';
% ft_databrowser(cfg)
% 
% % after ica
% cfg                         = [];
% cfg.channel                 = [{'megplanar'},strcat('-',subjectdata.badchannels)];
% cfg.continuous              = 'yes';
% cfg.blocksize               = 100;
% cfg.artfctdef.threshold.artifact = artifact_threshold1;
% cfg.plotevents              = 'no';
% ft_databrowser(cfg,data_preprocessed{1})
% 



