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
subjects = {'subject19'};

% choose files
files2preproc = 'stories_maxfilter';

% bandpass fequency
bp_freq = [4,30];

% apply ica
ica_on = 0;

% choose files for detected ica components
ica_files = 'all_maxfilter';

% downsample data
downsample_data = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
old = 'Z:'; new = '/media/till/Samsung_T5';
new = 'Z:';

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
    
    %% again trial rejection for the rest
    % better to do it after removing of the hugh ecg artifacts
    % magnetometer
%     cfg                              = [];
%     cfg.artfctdef.threshold.channel  = 'megmag'; 
%     %cfg.artfctdef.threshold.range    = 2500*10^-15; % 1000 fT (Stefan,Rupp)
%     cfg.artfctdef.threshold.min      = -1000*10^-15;
%     cfg.artfctdef.threshold.max      = 1000*10^-15;
%     cfg.artfctdef.threshold.bpfilter = 'no';
%     [~, artifact_threshold1]         = ft_artifact_threshold(cfg,data);
%     
%     % gradiometer
%     cfg                              = [];
%     cfg.artfctdef.threshold.channel  = 'megplanar'; 
%     %cfg.artfctdef.threshold.range    = 2500*10^-15/(4*10^-2); % 800 fT (Stefan,Rupp)
%     cfg.artfctdef.threshold.min      = -1000*10^-15/(4*10^-2);
%     cfg.artfctdef.threshold.max      = 1000*10^-15/(4*10^-2);
%     cfg.artfctdef.threshold.bpfilter = 'no';
%     [~, artifact_threshold2]         = ft_artifact_threshold(cfg,data);

% %     z-value, std
%     cfg                             = [];
%     cfg.artfctdef.zvalue.channel    = 'megplanar';
%     cfg.artfctdef.zvalue.cutoff     = 6;
%     cfg.artfctdef.zvalue.trlpadding = 0;
%     cfg.artfctdef.zvalue.artpadding = 0;
%     cfg.artfctdef.zvalue.fltpadding = 0;
%     cfg.artfctdef.zvalue.interactive = 'yes';
%     [cfg1,artifact] = ft_artifact_zvalue(cfg,data);
%     
%     cfg                              = [];
%     cfg.artfctdef.threshold.artifact = [artifact_threshold1;artifact_threshold2];
%     cfg.artfctdef.zvalue.artifact    = artifact;
%     data1 = ft_rejectartifact(cfg,data);
    
    %% baseline correction
    cfg                = [];
    cfg.baselinewindow = [-0.5 0];
    cfg.demean         = 'yes';
    data               = ft_preprocessing(cfg,data); 
    
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
data_preprocessed1     = ft_appenddata(cfg,data_preprocessed{:});
data_preprocessed1.hdr = hdr;
% 
% cfg         = [];
% cfg.method  = 'summary';
% cfg.channel = 'meg';
% data_clean  = ft_rejectvisual(cfg, data_preprocessed1);

cfg                             = [];
cfg.artfctdef.zvalue.channel    = 'all'; %['megmag',append('-',subjectdata.badchannels)];
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;
cfg.artfctdef.zvalue.interactive = 'yes';
cfg.artfctdef.zvalue.cutoff     = 5;
cfg.artfctdef.zvalue.abs = 'yes';
%cfg.artfctdef.zvalue.rectify       = 'yes';
[cfg1,artifact_high] = ft_artifact_zvalue(cfg,data_preprocessed1);

cfg                           = [];
cfg.artfctdef.zvalue.artifact = [artifact_high];
data_preprocessed1            = ft_rejectartifact(cfg,data_preprocessed1);







end

%% average data
cfg      = [];
data_avg = ft_timelockanalysis(cfg,data_preprocessed);

%% combine planar gradient
cfg            = [];
cfg.method     = 'sum';
cfg.updatesens = 'yes'; % need old information for dipole fitting 'no'
data_avg_cmb   = ft_combineplanar(cfg,data_avg);

% magnetometer
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag.lay';
cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,data_avg);
suptitle('avergage magnetometer')

figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar.lay';
cfg.xlim       = [-0.5, 0.5];
ft_multiplotER(cfg,data_avg);
suptitle('avergage gradiometer')

% planar gradient
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb.lay';
ft_multiplotER(cfg, data_avg_cmb);
suptitle('average combined gradiometers')






%% Clean up
rmpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new))
rmpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new))



cfg                              = [];
cfg.artfctdef.threshold.channel  = 'megmag'; 
%cfg.artfctdef.threshold.range    = 2000*10^-15; % 1000 fT (Stefan,Rupp)
cfg.artfctdef.threshold.min      = -1000*10^-15;
cfg.artfctdef.threshold.max      = 1000*10^-15;
cfg.artfctdef.threshold.bpfilter = 'no';
[~, artifacts]                   = ft_artifact_threshold(cfg,data_preprocessed);


cfg = [];
cfg.artfctdef.threshold.artifact = artifacts;
cfg.channel = 'megmag';
cfg.viewmode = 'butterfly';
cfg.ylim     = [-1000*10^-15,1000*10^-15];
ft_databrowser(cfg,data_preprocessed)




cfg                              = [];
cfg.artfctdef.threshold.channel  = 'megplanar'; 
%cfg.artfctdef.threshold.range    = 800*10^-15/(4*10^-2); % 1000 fT (Stefan,Rupp)
cfg.artfctdef.threshold.min      = -1000*10^-15/(4*10^-2);
cfg.artfctdef.threshold.max      = 1000*10^-15/(4*10^-2);
cfg.artfctdef.threshold.bpfilter = 'no';
[~, artifacts]                   = ft_artifact_threshold(cfg,data_preprocessed);

cfg = [];
cfg.artfctdef.threshold.artifact = artifacts;
cfg.channel = 'megplanar';
cfg.viewmode = 'butterfly';
cfg.ylim     = [-1000*10^-15/(4*10^-2),1000*10^-15/(4*10^-2)];
ft_databrowser(cfg,data_preprocessed)


