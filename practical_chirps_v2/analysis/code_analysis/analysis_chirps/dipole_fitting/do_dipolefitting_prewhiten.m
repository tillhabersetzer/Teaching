 close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to give you a first overview about the averages and
% dipolefitting
% does baseline normalization
%
% with prewhitening!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
% subjects = {'subject04'};
for i = 3:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjects{i-2} = [subject,num2str(i)]; 
end

% channels of interest for sourcemodelling
coi = 'meg';

% timewindow of interest for dipolefitting
timewin = [0.05 0.15];

% check results in plots
check = 0;

% choose files
% files2preproc = 'stories_maxfilter';

% bandpass fequency
% bp_freq = [1, 40];

% apply ica
ica_status = 1;

% choose files for detected ica components
% ica_files = 'all_maxfilter';

% apply baseline correction
% baseline_correction_status = 1;

% downsample data
% downsample_status = 0;
% fs_down           = 200;

% save data
save_data = 1;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);
% addpath for ica functions
addpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);
% addpath for preprocessing function
addpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);

if ica_status
    add = '_ica';
else 
    add = '';
end

%% analysis

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj

%% preprocessing
% config.subject                    = subjects{s};
% config.files2preproc              = files2preproc;
% config.bp_freq                    = bp_freq;
% config.ica_status                 = ica_status;
% config.ica_files                  = ica_files;
% config.baseline_correction_status = baseline_correction_status;
% config.downsample_status          = downsample_status;
% config.fs_down                    = fs_down;
% 
% data_preprocessed = preprocess_chirps(config);
eval(subjects{s})
data              = importdata(fullfile(subjectdata.chirps_sensorlevel,[subjectdata.subjectname,'_averages',add,'.mat']));
config            = data.info;
data_preprocessed = data.data_preprocessed;
clear data;
    
%% Noise-covariance estimation
% for a correct noise-covariance estimation it is important that 
% you used the cfg.demean = 'yes';

filenames_noise = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
noise           = give_noise(filenames_noise,subjectdata);

% if config.ica_status
%    noise = reject_independent_components(noise,subjectdata,config.ica_files); 
% end

cfg                  = [];
cfg.channel          = coi;
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg,noise);
avg                  = ft_timelockanalysis(cfg,data_preprocessed);

if check
    % noise covariance matrix
    %------------------------
    selmag  = ft_chantype(avg_noise.label, 'megmag');
    selgrad = ft_chantype(avg_noise.label, 'megplanar');
    C = avg_noise.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    figure
    subplot(1,2,1)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('MEG sensor covariance matrix')
    
    % singular values
    %----------------
    [~,s,~] = svd(avg_noise.cov);
    subplot(1,2,2)
    plot(log10(diag(s)),'o');
    grid on
    title('Singular values of a MEG sensor covariance matrix')
    sgtitle(subjectdata.subjectname)
end

%% prewhitening
% the following lines detect the location of the first large 'cliff' in the 
% singular value spectrum of the grads and mags
switch coi
    case 'meg'
        sensors = {'megmag','megplanar'};
    case {'megplanar','megmag'}
        sensors = {coi};
end
kappa_noise = give_kappa_value(avg_noise.cov,avg_noise.label,sensors);
kappa_data  = give_kappa_value(avg.cov,avg.label,sensors);

kappa_data(kappa_data<50) = [];
kappa                     = min([kappa_noise,kappa_data]); 
if kappa < 50; error('unexpected kappa value'); end

if ica_status % necessary for ft_denoise to work
    [i1,i2] = match_str(data_preprocessed.grad.label,avg_noise.grad.label);
    data_preprocessed.grad.chanunit(i1) = avg_noise.grad.chanunit(i2);
    data_preprocessed.grad.chantype(i1) = avg_noise.grad.chantype(i2);
end

cfg                 = [];
cfg.channel         = coi;
cfg.kappa           = kappa; % ensures use of regularized inverse
data_preprocessed_w = ft_denoise_prewhiten(cfg, data_preprocessed, avg_noise);
noise_w             = ft_denoise_prewhiten(cfg, noise, avg_noise);

cfg                  = [];
cfg.removemean       = 'yes'; % default for covariance computation
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_w                = ft_timelockanalysis(cfg,data_preprocessed_w);

% check covariance
%-----------------
if check 
    cfg                  = [];
    cfg.channel          = coi;
    cfg.covariance       = 'yes';
    cfg.covariancewindow = 'all';
    avg_noise_w          = ft_timelockanalysis(cfg,noise_w);

    % noise and data covariance matrix after whitening
    %-------------------------------------------------
    selmag  = ft_chantype(avg_noise_w.label, 'megmag');
    selgrad = ft_chantype(avg_noise_w.label, 'megplanar');
    C = avg_noise_w.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    figure
    subplot(2,2,1)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('noise')
    
    selmag  = ft_chantype(avg_w.label, 'megmag');
    selgrad = ft_chantype(avg_w.label, 'megplanar');
    C = avg_w.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    subplot(2,2,2)
    imagesc(C);
    hold on;
    plot(102.5.*[1 1],[0 306],'w','linewidth',2);
    plot([0 306],102.5.*[1 1],'w','linewidth',2);
    title('data')
    
    % singular values after whitening
    %--------------------------------
    [~,sv,~] = svd(avg_noise_w.cov);
    subplot(2,2,3)
    plot(log10(diag(sv)),'o');
    grid on
    title('noise')
    
    [~,sv,~] = svd(avg_w.cov);
    subplot(2,2,4)
    plot(log10(diag(sv)),'o');
    grid on
    title('data')
    sgtitle([subjectdata.subjectname,': after whitening'])
end

%% Fit a dipole model to the MEG data
%--------------------------------------------------------------------------
headmodel     = importdata(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_headmodel.mat'])); % mm
sourcemodel   = importdata(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_grid.mat'])); % mm
mri_segmented = importdata(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_mri_segmented.mat'])); % mm
template_grid = importdata(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ...
                            'template\sourcemodel\standard_sourcemodel3d10mm.mat']);
template_grid = ft_convert_units(template_grid,'mm');

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

% prepare leadfield
%------------------
cfg                       = [];
cfg.grad                  = data_preprocessed_w.grad;  % sensor information
cfg.channel               = coi;                       % the used channels
cfg.sourcemodel           = sourcemodel;               % source points
cfg.headmodel             = headmodel;                 % volume conduction model
cfg.singleshell.batchsize = 5000;                      % speeds up the computation
leadfield                 = ft_prepare_leadfield(cfg); % NOTE: input of the whitened data ensures the correct sensor definition to be used.

% symmetrical grid search
%------------------------
cfg             = [];
cfg.latency     = timewin;
cfg.numdipoles  = 2;
cfg.symmetry    = 'x';
cfg.sourcemodel = leadfield;
cfg.gridsearch  = 'yes';
cfg.headmodel   = headmodel;
cfg.senstype    = 'meg';
cfg.channel     = 'meg';
source          = ft_dipolefitting(cfg, avg_w);

% positions in mni space
idx         = dsearchn(sourcemodel.pos,source.dip.pos);
pos_sym_mni = template_grid.pos(idx,:);

% check
%------
if check

figure
hold on
ft_plot_dipole(source.dip.pos(1,:), mean(source.dip.mom(1:3,:),2), 'color', 'b', 'unit', 'mm')
ft_plot_dipole(source.dip.pos(2,:), mean(source.dip.mom(4:6,:),2), 'color', 'b', 'unit', 'mm')

pos = mean(source.dip.pos,1);
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
cfg.latency         = timewin;
cfg.numdipoles      = 2;
cfg.symmetry        = [];
cfg.gridsearch      = 'no';
cfg.dip.pos         = source.dip.pos;
cfg.sourcemodel     = leadfield;
cfg.headmodel       = headmodel;
cfg.channel         = 'meg';
cfg.senstype        = 'meg';
source_nosym = ft_dipolefitting(cfg,avg_w);

% positions in mni space
idx           = dsearchn(sourcemodel.pos,source_nosym.dip.pos);
pos_nosym_mni = template_grid.pos(idx,:);

% check
%------
if check 

figure
hold on

ft_plot_dipole(source.dip.pos(1,:), mean(source.dip.mom(1:3,:),2), 'color', 'g', 'unit', 'mm')
ft_plot_dipole(source.dip.pos(2,:), mean(source.dip.mom(4:6,:),2), 'color', 'g', 'unit', 'mm')

ft_plot_dipole(source_nosym.dip.pos(1,:), mean(source_nosym.dip.mom(1:3,:),2), 'color', 'm',  'unit', 'mm')
ft_plot_dipole(source_nosym.dip.pos(2,:), mean(source_nosym.dip.mom(4:6,:),2), 'color', 'm',  'unit', 'mm')

pos = mean(source.dip.pos,1);
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
cfg.sourcemodel = leadfield;
cfg.headmodel   = headmodel;
cfg.channel     = 'meg';
cfg.senstype    = 'meg';
source_all      = ft_dipolefitting(cfg, avg_w); % estimate the amplitude and orientation

% check
%------
if check 

figure
c = {'b','g','r'};
subplot(2,1,1); title('meg: probably right')
hold on
%plot(source_all.time, source_all.dip.mom(1:3,:), '-','Color',{'r','y','b'})
arrayfun(@(i) plot(source_all.time,source_all.dip.mom(i,:),'-','color',c{i}),1:3)
legend({'x', 'y', 'z'});
%axis([-0.1 0.4 -40e-3 40e-3])
grid on

subplot(2,1,2); title('meg: probably left')
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
% cfg             = [];
% cfg.model       = 'moving';  % default is rotating
% cfg.latency     = timewin; % [0.050 0.20]
% cfg.numdipoles  = 2;
% cfg.gridsearch  = 'no';
% cfg.dip.pos     = source_nosym.dip.pos;
% cfg.sourcemodel = leadfield;
% cfg.headmodel   = headmodel;
% cfg.channel     = 'meg';
% cfg.senstype    = 'meg';
% sources         = ft_dipolefitting(cfg, avg_w);


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

info                  = config;
info.prewhitening     = 'yes';
info.number_of_epochs = length(data_preprocessed.trial);
info.kappa_noise      = kappa_noise;
info.kappa_data       = kappa_data;
info.kappa            = kappa;
info.kappa_order      = sensors;

results.source            = source;       % fitted for specified latency
results.source_nosym      = source_nosym; % fitted for specified latency
results.source_timeseries = source_all;   % fitted for complete interval
% results.moving_source     = sources;    % moving dipole
results.pos_sym_mni       = pos_sym_mni;
results.pos_nosym_mni     = pos_nosym_mni;
save(fullfile(subjectdata.chirps_dipolefitting,[subjectdata.subjectname,'_results_dipolefitting',add,'.mat']),'results')

clear results

% template test
%--------------------------------------------------------------------------
% template_grid = importdata(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ...
%                             'template\sourcemodel\standard_sourcemodel3d10mm.mat']);
% template_grid = ft_convert_units(template_grid,'mm');
% template_mri  = ft_read_mri(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ... % mm
%                             'template\anatomy\single_subj_T1.nii']);
% atlas = ft_read_atlas('Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\template\atlas\aal\ROI_MNI_V4.nii');
% 
% idx     = dsearchn(sourcemodel.pos,source_nosym.dip.pos);
% pos_mni = template_grid.pos(idx,:);
% 
% cfg       = [];
% cfg.atlas = atlas;
% cfg.roi   = pos_mni;
% labels    = my_atlas_lookup(atlas,pos_mni(2,:),'coordsys','mni');
% 
% source_mni         = source_nosym;
% source_mni.dip.pos = pos_mni;
% 
% cfg               = [];
% cfg.method        = 'ortho';
% % cfg.funparameter  = 'mom';
% %cfg.location      = [64 -32 8];
% %cfg.funcolormap   = 'jet';
% %cfg.latency       = 0.1
% %cfg.avgovertime   = 'yes';
% ft_sourceplot(cfg,source_mni);
%     
% figure
% hold on
% 
% ft_plot_dipole(source_mni.dip.pos(1,:), mean(source_mni.dip.mom(1:3,:),2), 'color', 'g', 'unit', 'mm')
% ft_plot_dipole(source_mni.dip.pos(2,:), mean(source_mni.dip.mom(4:6,:),2), 'color', 'g', 'unit', 'mm')
% 
% pos = mean(source_mni.dip.pos,1);
% ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
% ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
% ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)
% 
% ft_plot_crosshair(pos, 'color', [1 1 1]/2);
% axis tight
% axis off
% view(12, -10)
% 
% % source.dip.pos
% % idx = dsearchn(sourcemodel.pos,source.dip.pos)
% 
% sourcemodel.pos(idx(2),:)
%--------------------------------------------------------------------------

end % subjects

disp('do_dipolefitting finished')

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files']);
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);
rmpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);