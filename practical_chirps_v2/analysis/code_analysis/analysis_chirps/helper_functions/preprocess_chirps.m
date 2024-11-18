function data_preprocessed = preprocess_chirps(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing of chirp data
% need to add following paths manually
% 
% addpath for subject_files information
% addpath(replace(['Z:' filesep 'analysis' filesep 'subject_files'],old,new));
% addpath for ica functions
% addpath(replace(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'],old,new));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% subject 
subject = config.subject;
% file
files2preproc = config.files2preproc;
% bandpass fequency
bp_freq = config.bp_freq;
% apply ica
ica_status = config.ica_status;
% choose files for detected ica components
ica_files = config.ica_files;
% apply baseline correction
baseline_correction_status = config.baseline_correction_status;
% downsample data
downsample_status = config.downsample_status;
fs_down           = config.fs_down;
%--------------------------------------------------------------------------

%% preprocessing
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
    if ica_status 
        data = reject_independent_components(data,subjectdata,ica_files);
    end
    
    %% filter data
    cfg              = [];
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = bp_freq;
    cfg.dftfilter    = 'yes';        % enable notch filtering to eliminate power line noise
    cfg.dftfreq      = [50 100 150]; % set up the frequencies for notch filtering
    cfg.dftreplace   = 'neighbour';
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
    cfg                  = [];
    cfg.trl              = trl_new;          
    data_preprocessed{f} = ft_redefinetrial(cfg,data); 
    
    %% again trial rejection for the rest
    % magnetometer
%     cfg                              = [];
%     cfg.artfctdef.threshold.channel  = 'megmag'; 
%     cfg.artfctdef.threshold.min      = -1000*10^-15;
%     cfg.artfctdef.threshold.max      = 1000*10^-15;
%     cfg.artfctdef.threshold.bpfilter = 'no';
%     [~, artifact_threshold1]         = ft_artifact_threshold(cfg,data);
%     
%     % gradiometer
%     cfg                              = [];
%     cfg.artfctdef.threshold.channel  = 'megplanar'; 
%     cfg.artfctdef.threshold.min      = -1000*10^-15/(4*10^-2);
%     cfg.artfctdef.threshold.max      = 1000*10^-15/(4*10^-2);
%     cfg.artfctdef.threshold.bpfilter = 'no';
%     [~, artifact_threshold2]         = ft_artifact_threshold(cfg,data);
%     
%     cfg                              = [];
%     cfg.artfctdef.threshold.artifact = [artifact_threshold1;artifact_threshold2];
%     data = ft_rejectartifact(cfg,data);
    
end % files
    
%% append data
hdr                   = data_preprocessed{1}.hdr;
cfg                   = [];
cfg.keepsampleinfo    = 'no';
data_preprocessed     = ft_appenddata(cfg,data_preprocessed{:});
data_preprocessed.hdr = hdr;

%% reject bad epochs based on z-value
cfg                              = [];
cfg.artfctdef.zvalue.channel     = 'all'; %['megmag',append('-',subjectdata.badchannels)];
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0;
cfg.artfctdef.zvalue.fltpadding  = 0;
cfg.artfctdef.zvalue.interactive = 'no';
cfg.artfctdef.zvalue.cutoff      = 5;
[~,artifacts]                    = ft_artifact_zvalue(cfg,data_preprocessed);

cfg                           = [];
cfg.artfctdef.zvalue.artifact = artifacts;
data_preprocessed             = ft_rejectartifact(cfg,data_preprocessed);

%% baseline correction
if baseline_correction_status
    cfg                = [];
    cfg.baselinewindow = [-0.5 0];
    cfg.demean         = 'yes';
    data_preprocessed  = ft_preprocessing(cfg,data_preprocessed);
end

%% downsample data
if downsample_status
    cfg               = [];
    cfg.resamplefs    = fs_down;
    cfg.detrend       = 'no';
    data_preprocessed = ft_resampledata(cfg,data_preprocessed);
end

end