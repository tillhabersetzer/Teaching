close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs identification of artifact-related independent
% components (part 4)
% - performs coherence analysis to identify independent components
% - maybe correlation? crosscorelation?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subjects
subjectlist = {'subject24'};

% choose files
% files2preproc = {'stories_maxfilter','olsa_maxfilter'};
files2preproc = {'all_maxfilter'};
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
% data containing artifact trials and independent components of sensor data
filename_data           = [subjectdata.subjectname '_' files2preproc{m} '_artifact_component_data_4ica.mat'];
artifact_component_data = importdata([path2save filesep filename_data]);
             
num_artifact   = length(artifact_component_data.arti_comp_data);
label_artifact = artifact_component_data.label;

%% 2.) check for artefact-related independent components - coherence
%--------------------------------------------------------------------------
%% 2.1) compute a frequency decomposition of trials for all artifacts and sensors
freq       = cell(1,num_artifact);
cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [0 50];
cfg.taper  = 'hanning';
cfg.pad    = 'maxperlen';

for a = 1: num_artifact      
       freq{a} = ft_freqanalysis(cfg,artifact_component_data.arti_comp_data{a});      
end

%% 2.2) compute coherence between all independent components and the artifacts
fdcomp         = cell(1,num_artifact);
cfg            = [];
cfg.jackknife  = 'no';
cfg.method     = 'coh';

for a = 1:num_artifact 
    cfg.channelcmb = {label_artifact{a},'all'};   
    fdcomp{a}      = ft_connectivityanalysis(cfg, freq{a});  
end

%% 2.3) look at the coherence spectrum between all components and the artifacts

% channels of interest for artefact detection
coi_artifact = {'ECG','EOG001','EOG002'};
num_artifact = length(coi_artifact);

% choose here - values correspond to index in coi_artifact/coi_sensor
%--------------------------------------------------------------------------
artifact = [1,2,3]; % ECG etc.
%--------------------------------------------------------------------------
% find index of elements in arti_comp_data to loop over
ind_a  = contains(label_artifact,coi_artifact(artifact),'IgnoreCase',true)';
a_loop = find(ind_a);

for a = a_loop
    fig_name = [subjectlist{s} ' | ' files2preproc{m} ' | Sensor: meg | Artifact: ',label_artifact{a}];
    figure('name',fig_name);
    subplot(2,1,1); plot(fdcomp{a}.freq, abs(fdcomp{a}.cohspctrm));
    subplot(2,1,2); imagesc(abs(fdcomp{a}.cohspctrm));
    suptitle(fig_name)
end  

%% 3.) try pearson correlation
%--------------------------------------------------------------------------
% fieldtrip stuff does not work
% -> apply simple correlation coefficient, therefore you assume no shift
%    between the components and the artifacts

%% 3.1) compute correlations
pearson_corr = cell(1,num_artifact);

for a = 1: num_artifact 
    % cfg.channelcmb = {label_artifact{a},'all' ,};
    % number of trials
    num_trials = length(artifact_component_data.arti_comp_data{a}.trial);
    % number of components
    num_comp   = length(artifact_component_data.arti_comp_data{a}.label)-1;
    % initilaize matrix for correlations
    correlations = zeros(num_comp,num_trials);
    for i = 1:num_trials
        % data for correlattion
        data2corr         = artifact_component_data.arti_comp_data{a}.trial{i};
        correlations(:,i) = corr(data2corr(1,:)',data2corr(2:end,:)','Type','Pearson')'; % columswise comparison
    end
    pearson_corr{a} = correlations;
    clear num_trials correlations
end      

%% 3.2) plot correlations

% channels of interest for artefact detection
coi_artifact = {'ECG','EOG001','EOG002'};

% choose here - values correspond to index in coi_artifact/coi_sensor
%--------------------------------------------------------------------------
artifact = [1,2,3]; % ECG etc.
%--------------------------------------------------------------------------
% find index of elements in arti_comp_data to loop over
ind_a  = contains(label_artifact,coi_artifact(artifact),'IgnoreCase',true)';
a_loop = find(ind_a);

for a = a_loop
    fig_name = [subjectlist{s} ' | ' files2preproc{m} ' | Sensor: meg | Artifact: ',label_artifact{a}];
    figure('name',fig_name);
    % average correlation coefficient over all trials for one
    % component
    subplot(2,1,1); plot(mean(pearson_corr{a},2),'-o');
    xlabel('component')
    ylabel('mean correlation coefficient')
    title('mean correlation coefficient for component')
    subplot(2,1,2); imagesc(pearson_corr{a});
    xlabel('trials')
    ylabel('components')
    suptitle(fig_name)
end  

         
%% 4.) check for artifact-related independent components - timelocked/averaged components
%--------------------------------------------------------------------------

%% 4.1) average the components timelocked to the artifact-complex
% only valid for ecg artifacts cause the were gathered timelocked and all
% trials have the same length
cfg      = [];
ind_a    = contains(label_artifact,{'ECG'},'IgnoreCase',true);
a        = find(ind_a);
 
timelock = ft_timelockanalysis(cfg,artifact_component_data.arti_comp_data{a});
    
%% 4.2) look at the timelocked/averaged components and compare them with the ECG

fig_name = [subjectlist{s} ' | ' files2preproc{m} ' | Sensor: meg | Artifact: ',label_artifact{a}];
figure('name',fig_name);
title('comparison of timelocked components and ecg')
subplot(2,1,1); plot(timelock.time, timelock.avg(1,:))
subplot(2,1,2); plot(timelock.time, timelock.avg(2:end,:))
suptitle(fig_name)
figure('name',fig_name);
title('comparison of timelocked components and ecg')
subplot(2,1,1); plot(timelock.time, timelock.avg(1,:))
subplot(2,1,2); imagesc(timelock.avg(2:end,:));
suptitle(fig_name)

% asks user for input
response = input('Next file? [y]/[n]:','s');
if strcmpi(response,'n')
    break
end

end % files

% aks user for input
response = input('Next subject? [y]/[n]:','s');
if strcmpi(response,'n')
    break
end

end % subjects

rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])



