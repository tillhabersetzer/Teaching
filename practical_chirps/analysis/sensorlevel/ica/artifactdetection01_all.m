close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs identification of artifact-related independent
% components (part 1)
% Therefore it is possible to use multiple files to compute independent 
% component analysis on
% - it calculates the decomposition of meg data into it's 
%   independent components
% - it also saves continous preprocessed data which is needed for the later
%   identification of artifact-related independent components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subjects
% subjectlist = {'subject02'};
for i = 16:24 %24
    if i<10; subject='subject0'; else subject='subject'; end
    subjectlist{i-15} = [subject,num2str(i)]; 
end


% choose files
% files2preproc = {'stories_maxfilter','olsa_maxfilter'};
files2preproc = {'all_maxfilter'};

% option to reject trials before ica
trialrejection = 1;

% option to save 
save_data = 0;

% option to clear data
c = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);
% addpath for ica functions
addpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions']);
% addpath for preprocessing function
addpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);

%% initialize

% loop over subjects 
%-------------------
N_subj  = length(subjectlist);
N_files = length(files2preproc);

for s = 1:N_subj
    
% subject selected out of codelist 
subject = subjectlist{s};
eval(subject)

% loop over selected files
%-------------------------

for m = 1:N_files % if you want to separate between olsa and stories

filenames = get_filenames(subjectdata,files2preproc{m},1); % 1: remove badfiles
F         = length(filenames);

% these cells are saved for further processing
data_sensor   = cell(1,F);
data_artifact = cell(1,F);
triallabel    = zeros(F,2); % filenames x [start,end]
% needed for ica on all files
data_epoched  = cell(1,F);


%% check rank of data
[ranks,ranks_table] = compute_ranks1(subjectdata,filenames); 
r                   = min(ranks,[],'all'); 
%% preprocessing
for f = 1:F % loop over files
   
    %% 1.) load and filter continuous data
    cfg              = [];
    cfg.dataset      = [subjectdata.rawdatadir filesep filenames{f} '.fif'];
    % channels of interest and exclude bad channels
    % cfg.channel      = [{'meg','eog','ecg'},strcat('-',subjectdata.badchannels)]; 
    cfg.channel      = {'meg','eog','ecg'}; 
    cfg.demean       = 'yes';
    cfg.detrend      = 'yes';
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;            % ensure that sensors are expressed in SI units
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = [0.1,150];
    cfg.dftfilter    = 'yes';        % enable notch filtering to eliminate power line noise
    cfg.dftfreq      = [50 100 150]; % set up the frequencies for notch filtering
    cfg.dftreplace   = 'neighbour';
    data_filtered    = ft_preprocessing(cfg);      

    %% 2.) separate the data
    coi_artifact     = {'ecg','eog'};
    coi_sensor       = 'meg';
    cfg              = [];
    cfg.channel      = coi_artifact;
    data_artifact{f} = ft_selectdata(cfg, data_filtered); 
    cfg.channel      = coi_sensor;
    data_sensor{f}   = ft_selectdata(cfg, data_filtered); 
    
    if c
        clear data_filtered 
    end
    
    %% 3.) define segment for ica
    cfg                      = [];
    cfg.dataset              = [subjectdata.rawdatadir filesep filenames{f} '.fif'];
    cfg.trialfun             = 'my_trialfun_ica';
    % select interval for ica dependet on trigger channels
    cfg.trialdef.eventtype   = 'STI101';
    cfg.trialdef.triallength = 5; 
    cfg                      = ft_definetrial(cfg);
    trl                      = cfg.trl;

    %% 4.) epoch data
    cfg             = [];
    % multiple segmentes between triggers
    cfg.trl         = trl;
    cfg.continuous  = 'yes';
    data_epoched{f} = ft_redefinetrial(cfg,data_sensor{f});
    
    % filter epoched data for line noise - only works for shorter segments
%     cfg             = [];
%     cfg.channel     = 'meg'; 
%     cfg.dftfilter   = 'yes';        % enable notch filtering to eliminate power line noise
%     cfg.dftfreq     = [50 100 150]; % set up the frequencies for notch filtering
%     data_epoched{f} = ft_preprocessing(cfg,data_epoched{f});
    
    cfg                      = [];
    % complete segments between triggers
    cfg.trl          = [trl(1,1),trl(end,2),0];
    cfg.continuous   = 'yes';
    data_artifact{f} = ft_redefinetrial(cfg,data_artifact{f});
    
    cfg            = [];
    % complete segments between triggers
    cfg.trl        = [trl(1,1),trl(end,2),0];
    cfg.continuous = 'yes';
    data_sensor{f} = ft_redefinetrial(cfg,data_sensor{f});
    
    %% 4.) downsample data to reduce memory load
    cfg              = [];
    cfg.resamplefs   = 400;
    cfg.detrend      = 'no';
    data_epoched{f}  = ft_resampledata(cfg,data_epoched{f}); 
    data_artifact{f} = ft_resampledata(cfg,data_artifact{f});
    data_sensor{f}   = ft_resampledata(cfg,data_sensor{f});
    
    %% 4.) clean ica data - reject heavily contaminated epochs    
    if trialrejection
        cfg                        = [];
        cfg.continuous             = 'yes';
        cfg.artfctdef.jump.channel = 'meg'; 
        [~, artifact_jump]         = ft_artifact_jump(cfg, data_epoched{f});

        cfg                        = [];
        cfg.continuous             = 'yes';
        cfg.artfctdef.clip.channel = 'all'; 
        [~, artifact_clip]         = ft_artifact_clip(cfg, data_epoched{f});

        cfg = [];
        % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
        cfg.artfctdef.reject        = 'partial';
        cfg.artfctdef.jump.artifact = artifact_jump;
        cfg.artfctdef.clip.artifact = artifact_clip;
        data_epoched{f}             = ft_rejectartifact(cfg,data_epoched{f});   
    end

end

%% 5.) prepare sensor data and artifact data for saving
hdr_sensor   = data_sensor{1}.hdr;
hdr_artifact = data_artifact{1}.hdr;

cfg                = [];
cfg.keepsampleinfo = 'no';
data_sensor        = ft_appenddata(cfg,data_sensor{:});
data_artifact      = ft_appenddata(cfg,data_artifact{:});

% add additional information to data
data_sensor.hdr          = hdr_sensor;
data_artifact.hdr        = hdr_artifact;
data_sensor.triallabel   = filenams;
data_artifact.triallabel = filenames;

if contains(files2preproc{m},'stories')
    path2save = subjectdata.ica_stories;
elseif contains(files2preproc{m},'olsa')
    path2save = subjectdata.ica_olsa;
elseif contains(files2preproc{m},'all')
    path2save = subjectdata.ica_all;
else
    error('Specified path2save does not exist!')
end 

% save artifact and sensor data
if save_data
    filename_new = [subjectdata.subjectname '_' files2preproc{m} '_data_artifact_4ica.mat'];
    save([path2save filesep filename_new],'data_artifact','-v7.3')
    filename_new = [subjectdata.subjectname '_' files2preproc{m} '_data_sensor_4ica.mat'];
    save([path2save filesep filename_new],'data_sensor','-v7.3')   
end

if c
    clear data_artifact data_sensor 
end

%% 6.) concatenate data of all files for ica
cfg                = [];
cfg.keepsampleinfo = 'no';
data_epoched       = ft_appenddata(cfg,data_epoched{:});

%% 7.) ica on different sensortypes
% downsample data for ICA 
% include sensor channel, not artifact channels
% cfg              = [];
% cfg.resamplefs   = 200;
% cfg.detrend      = 'no';
% data_epoched = ft_resampledata(cfg,data_epoched);       

% perform ica
cfg                 = [];
% perform pca for maxfiltered meg data
%cfg.runica.pca      = 64; % see natmeg tutorial for combined meg/eeg
cfg.runica.pca      = r; % rank of data
cfg.runica.sphering = 'on'; % is default (for different sensortypes, variances are equalized)
cfg.method          = 'runica';
cfg.channel         = 'meg';
components          = ft_componentanalysis(cfg,data_epoched);  

%% 8.) save ica data and ranks

if save_data
    filename_new = [subjectdata.subjectname '_' files2preproc{m} '_ica_comp.mat'];
    save([path2save filesep filename_new],'components','-v7.3')

    filename_new = [subjectdata.subjectname '_' files2preproc{m} '_ranks.mat'];
    save([path2save filesep filename_new],'ranks_table','-v7.3')
end

end % files
end % subjects

% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])


% only possible when data has not been cleared
%--------------------------------------------------------------------------
% % show also recorded artifact channels
% cfg                         = [];
% cfg.channel                 = {'eog','ecg'}; % components to be plotted
% cfg.blocksize               = 750;
% cfg.plotevents              = 'no';
% cfg.viewmode                = 'vertical';
% ft_databrowser(cfg,data_artifact);
% 
% data_artifact.hdr = data_filtered.hdr;
% % show also recorded artifact channels
% cfg                         = [];
% %cfg.channel                 = {'eog','ecg'}; % components to be plotted
% cfg.channel                 = 'megmag'; % components to be plotted
% cfg.blocksize               = 750;
% cfg.plotevents              = 'no';
% ft_databrowser(cfg,data_sensor);



%% functions
%--------------------------------------------------------------------------
% both functions seems to give the same results!

function [ranks,ranks_table] = compute_ranks1(subjectdata,filenames)
sensors = {'megmag','megplanar'};
ranks   = zeros(length(filenames),length(sensors));

for f = 1:length(filenames)   
cfg                  = [];
cfg.dataset          = [subjectdata.rawdatadir filesep filenames{f} '.fif'];
cfg.coilaccuracy     = 0; 
data                 = ft_preprocessing(cfg);
 
cfg                  = [];
cfg.channel          = 'meg';
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg                  = ft_timelockanalysis(cfg,data);
kappa                = give_kappa_value(avg.cov,avg.label,sensors);       
ranks(f,:)           = kappa;
end

ranks_table = array2table(ranks,'RowNames',filenames,'VariableNames',sensors);
   
end

function [ranks,ranks_table] = compute_ranks2(subjectdata,filenames)
% computes the rank of the matrix containing the specified sensor-channels
% here you can see that maxfiltered data is rank-deficient 
sensors = {'megmag','megplanar','meg'};
ranks   = zeros(length(filenames),length(sensors));
for f = 1:length(filenames)   
cfg              = [];
cfg.dataset      = [subjectdata.rawdatadir filesep filenames{f} '.fif'];
cfg.coilaccuracy = 0; 
data             = ft_preprocessing(cfg);
    for s = 1:length(sensors)
        cfg         = [];
        %cfg.channel = cat(2,sensors{s},strcat('-',subjectdata.badchannels));
        cfg.channel = sensors{s};
        data_select = ft_selectdata(cfg,data);
        ranks(f,s)  = rank(data_select.trial{1}*data_select.trial{1}');
    end
end

ranks_table = array2table(ranks,'RowNames',filenames,'VariableNames',sensors);

end


