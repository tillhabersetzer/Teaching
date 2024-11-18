close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs identification of artifact-related independent
% components (part 3)
% - seperates sensor channel data (EEG/MEG) into artifact related trials
% - then it decomposes the trials into it's independent components based 
%   on already computed unmixing matrices (artifactdetection01) to relate
%   this independent components to the artifact channels.
% - therefore it generates a matrix containing independent component trials
%   and artifact trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subjects
% subjectlist = {'subject04'};
for i = 16:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjectlist{i-15} = [subject,num2str(i)]; 
end


% choose files
%files2preproc = {'stories_maxfilter','olsa_maxfilter'};
files2preproc = {'all_maxfilter'};


% option to save 
save_data = 1;

% option to clear data
c = 1;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);

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

for m = 1:N_files

if contains(files2preproc{m},'stories')
    path2save = subjectdata.ica_stories;
elseif contains(files2preproc{m},'olsa')
    path2save = subjectdata.ica_olsa;
elseif contains(files2preproc{m},'all')
    path2save = subjectdata.ica_all;
else
    error('Specified path2save does not exist!')
end    

%% 1.) load data
% continuous data (meg,eeg)
filename_data = [subjectdata.subjectname '_' files2preproc{m} '_data_sensor_4ica.mat'];
data_sensor   = importdata([path2save filesep filename_data]);

% continous artifact data (eog,ecg)
filename_artifact = [subjectdata.subjectname '_' files2preproc{m} '_data_artifact_4ica.mat'];
data_artifact     = importdata([path2save filesep filename_artifact]);

% ica unmixing matrices
filename_ica = [subjectdata.subjectname '_' files2preproc{m} '_ica_comp.mat'];
components   = importdata([path2save filesep filename_ica]);
 
% artifact trials
filename_trials  = [subjectdata.subjectname '_' files2preproc{m} '_artifact_trials_4ica.mat'];
trials           = importdata([path2save filesep filename_trials]);

%% 2.) Separate channel by beforehand calculated artifact trials
coi_artifact = data_artifact.label;
num_artifact = length(coi_artifact);
% artifact channels
data_artifact_trl = cell(1,num_artifact);
% sensor channels
data_sensor_trl   = cell(1,num_artifact);

for a = 1: num_artifact
    
    % only channel of interest
    cfg                  = [];
    cfg.channel          = coi_artifact{a}; 
    data_artifact_trl{a} = ft_selectdata(cfg, data_artifact); 
    
    ind_trl              = contains(trials.label,coi_artifact{a});
    artifact_trials      = trials.artifact_trials{ind_trl};
    cfg                  = [];
    cfg.trl              = [artifact_trials, zeros(size(artifact_trials,1),1)];
    data_artifact_trl{a} = ft_redefinetrial(cfg,data_artifact_trl{a});
         
    % only data segment around artifact trials
    cfg                = [];
    cfg.trl            = [artifact_trials, zeros(size(artifact_trials,1),1)];
    data_sensor_trl{a} = ft_redefinetrial(cfg,data_sensor);

end

if c
    clear data_sensor data_artifact trials
end

%% 3.) resample to speed up the decomposition and frequency analysis, especially usefull for 1200Hz MEG data
% -> already resampled in artifactdetection01 to 400 Hz
% cfg            = [];
% cfg.resamplefs = 500;
% cfg.detrend    = 'no';  
% 
% for a = 1:num_artifact 
%     data_artifact_trl{a} = ft_resampledata(cfg,data_artifact_trl{a});       
%     data_sensor_trl{a}   = ft_resampledata(cfg,data_sensor_trl{a});       
% end

%% 4.) decompose the artifact-locked datasegments into components, using the previously found (un)mixing matrix
data_sensor_comp = cell(1,num_artifact);

cfg           = [];
cfg.unmixing  = components.unmixing;
cfg.topolabel = components.topolabel;
for a = 1: num_artifact
  data_sensor_comp{a} = ft_componentanalysis(cfg, data_sensor_trl{a});
end  

if c
    clear data_sensor_trl
end

%% 5.) append the artifact-channels to the sensor data
cfg                = [];
cfg.keepsampleinfo = 'no'; % because data originates from different datafiles

for a = 1:num_artifact
      data_sensor_comp{a} = ft_appenddata(cfg, data_artifact_trl{a},data_sensor_comp{a});
end

artifact_component_data.arti_comp_data = data_sensor_comp; % artifact x sensor cell
artifact_component_data.label          = coi_artifact;
artifact_component_data.overview       = cell2table(cell(size(data_sensor_comp)),'RowNames',{'meg'},...
                                                                     'VariableNames',coi_artifact);
if c
    clear data_artifact_trl data_sensor_comp
end

if save_data
    filename_new = [subjectdata.subjectname '_' files2preproc{m} '_artifact_component_data_4ica.mat'];
    save([path2save filesep filename_new],'artifact_component_data','-v7.3')
end

end % files

end % subjects

rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])

