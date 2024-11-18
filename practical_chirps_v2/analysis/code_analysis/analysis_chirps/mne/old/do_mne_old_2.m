close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!SCRIPT IS DESIGNED FOR LINUX!!!!
% to work on windows replace new = '/media/till/Samsung_T5' with new = old
% then everthing with the pays stays the same
%
% This script tries to give you a first overview about the averages and
% minimum norm estimates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject05'};

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [4, 30];

% apply ica
ica_on = 0;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% downsample data
downsample_data = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
old = 'Z:'; new = '/media/till/Samsung_T5';
%new = 'Z:';

% addpath for subject_files information
addpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new));
% addpath for ica functions
addpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new));

%% initialize

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj
    % subject selected
subject = subjects{s};
eval(subject)
subjectdata.rawdatadir = replace(subjectdata.rawdatadir,old,new);

%filenames = get_filenames(subjectdata,files2preproc);
filenames = get_filenames(subjectdata,files2preproc);

N_files   = length(filenames);

% loop over selected files
%-------------------------
data_preprocessed = cell(1,N_files);

%parfor f = 1:N_files
for f = 1:N_files
    
    path_dataset = replace([subjectdata.rawdatadir filesep filenames{f} '.fif'],old,new);
    % filter continuous data to avoid edge artifacts
    cfg              = [];
    cfg.dataset      = path_dataset;
    cfg.channel      = 'meg'; 
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;            % ensure that sensors are expressed in SI units
    data             = ft_preprocessing(cfg);   
    
    %% reject earlier specified independet components
    if ica_on  
        % change necessary paths
        subjectdata.ica_olsa    = replace(subjectdata.ica_olsa,old,new);
        subjectdata.ica_stories = replace(subjectdata.ica_stories, old, new);
        subjectdata.ica_all     = replace(subjectdata.ica_all,old,new);
 
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

%% Averaging and noise-covariance estimation
% for a correct noise-covariance estimation it is important that 
% you used the cfg.demean = ‘yes’

coi                  = 'megmag';

cfg                  = [];
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = [-inf 0]; % timepoints before the zero timepoints in trials
avg                  = ft_timelockanalysis(cfg, data_preprocessed);


% estimate noise covariance different (empty rooms)
empty_room = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
N          = length(empty_room);
noise      = cell(1,N);

for n = 1:N
    cfg              = [];
    cfg.dataset      = replace([subjectdata.rawdatadir filesep empty_room{n} '.fif'],old,new);
    cfg.channel      = coi; 
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;            % ensure that sensors are expressed in SI units
    noise{n}         = ft_preprocessing(cfg);   
end

hdr                = noise{1}.hdr;
cfg                = [];
cfg.keepsampleinfo = 'no';
noise              = ft_appenddata(cfg,noise{:});
noise.hdr          = hdr;

cfg                  = [];
cfg.channel          = coi;
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg,noise);

figure
subplot(1,2,1)
imagesc(avg_noise.cov);
title('empty room measurements')
subplot(1,2,2)
imagesc(avg.cov);
title('baseline interval')

%% Forward solution

% load headmodel
headmodel = importdata(replace(fullfile(subjectdata.headmodel,[subjectdata.subjectname,'_headmodel.mat']),old,new)); % mm
% load sourcemodel
sourcemodel = importdata(replace(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_sourcemodel_4k.mat']),old,new)); % mm

% prepare leadfield
%------------------
cfg                       = [];
cfg.grad                  = data_preprocessed.grad; % sensor information
cfg.channel               = coi;                    % the used channels
cfg.sourcemodel           = sourcemodel;            % source points
cfg.headmodel             = headmodel;              % volume conduction model
cfg.singleshell.batchsize = 5000;                   % speeds up the computation
leadfield                 = ft_prepare_leadfield(cfg);

% constrained minimum norm estimate to surface normal
%----------------------------------------------------
% nrm = normals(sourcemodel.pos, sourcemodel.tri,'vertex');
% % modify leadfield - single column leadfield
% for i = 1:length(nrm)
%     if ~isempty(leadfield.leadfield{i})
%         leadfield.leadfield{i} = leadfield.leadfield{i}*nrm(i,:)'; 
%     end
% end

% check
%------
% same units for plotting
dataset = replace([subjectdata.rawdatadir filesep 'fixer_olsa_tsss_mc.fif'],old,new);
shape   = ft_read_headshape(dataset,'unit','mm');
grad    = ft_convert_units(ft_read_sens(dataset,'senstype','meg'),'mm'); 

figure
ft_plot_headmodel(headmodel,'facealpha',0.5);
hold on
ft_plot_sens(grad, 'style', '*b');
ft_plot_headshape(shape);
ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', ...
             'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);

%% Inverse solution
% modify noise covariance estimation
avg.cov = avg_noise.cov;

cfg                    = [];
cfg.method             = 'mne';
cfg.sourcemodel        = leadfield;
cfg.headmodel          = headmodel;
cfg.mne.prewhiten      = 'yes';
cfg.mne.lambda         = 3; % scaling factor for noise
cfg.mne.scalesourcecov = 'yes';
cfg.mne.keepfilter     = 'yes';
source                 = ft_sourceanalysis(cfg,avg);

%% Visualize
signal = source.avg.pow(:,640); % 100ms
%signal = sqrt(source.avg.pow(:,640)); % moments are sqrt's of power estimates

figure
ft_plot_mesh(source, 'vertexcolor', signal);
view([180 0]);
h = light; 
set(h, 'position', [0 1 0.2]); 
lighting gouraud;
material dull;

% movie
cfg            = [];
%cfg.projectmom = 'yes';
sd             = ft_sourcedescriptives(cfg,source);

cfg              = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,sd);

% cfg                = [];
% cfg.latency        = 0.113;
% cfg.method         = 'surface';
% cfg.funparameter   = 'pow';
% cfg.maskparameter  = cfg.funparameter;
% %cfg.funcolorlim    = [0.0 1.2];
% cfg.funcolormap    = 'jet';
% %cfg.opacitylim     = [0.0 1.2];
% cfg.opacitymap     = 'rampup';
% cfg.projmethod     = 'nearest';
% cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
% %cfg.surfdownsample = 10;  % downsample to speed up processing
% ft_sourceplot(cfg, source);
% view ([90 0])             % rotate the object in the view


end

%% Clean up
rmpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new))
rmpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new))

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,avg);
suptitle('average magnetometer')

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
%cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,avg);
suptitle('average gradiometer')

% atlas information
%roi = ft_read_headshape('/media/till/Volume/mne_workbench/template_files/R.atlasroi.4k_fs_LR.shape.gii');
surf = ft_read_headshape('/media/till/Volume/mne_workbench/template_files/R.sphere.4k_fs_LR.surf.gii');

%label = ft_read_headshape('/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.L.aparc.4k_fs_LR.label.gii');
label1 = ft_read_atlas('/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.R.aparc.4k_fs_LR.label.gii');
label2 = ft_read_atlas('/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.R.aparc.a2009s.4k_fs_LR.label.gii');
%label = ft_read_atlas('/media/till/Volume/mne_workbench/template_files/R.atlasroi.4k_fs_LR.shape.gii');

% a = '/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.R.aparc.a2009s.4k_fs_LR.label.gii';
% b = '/media/till/Volume/freesurfer_subject_directory/er04er24_vp/workbench/er04er24_vp.R.aparc.4k_fs_LR.label.gii';
% label_a = ft_read_atlas({a, strrep(a,'.R.','.L.')});
% label_b = ft_read_atlas({b, strrep(b,'.R.','.L.')});







