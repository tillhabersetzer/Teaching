close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs identification of artifact-related independent
% components (part 2)
% - performs trial detection of artifacts in the recorded artifact channels
%   and saves it 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subjects
subjectlist = {'subject24'};

% choose files
% files2preproc = {'stories_maxfilter','olsa_maxfilter'};
files2preproc = {'all_maxfilter'};

% option to save
save_data = 1;

% option to plot found artifacts
check_artifacts = 1;

% choose additional artifact suppression if bipolar channels are
% contaminated by high valued artifacts
additional_suppression = 0;
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

%% find artifact trials
%--------------------------------------------------------------------------

% load data (for  detect_artifact_trials1)
%--------------------------------------------------------------------------
% continous artifact channels
filename_data = [subjectdata.subjectname '_' files2preproc{m} '_data_artifact_4ica.mat'];
data_artifact = importdata([path2save filesep filename_data]);
%--------------------------------------------------------------------------
num_artifact = length(data_artifact.label);
artifact_trl = cell(1,num_artifact);

if additional_suppression
% apply own artifact detection for curious events which spreads in all
% recorded bipolar channels
%--------------------------------------------------------------------------
    % peaks are biggest in ecg channels
    cfg                = [];
    cfg.threshold      = 1000; %1000;     % arbitrarily set threshold
    cfg.channel        = 'ECG003';
    cfg.pretim         = 1;        % sec / interval for artifact
    cfg.psttim         = 1;
    artifact_intervals = detect_threshold_artifacts(cfg,data_artifact);

    % clean all artifact channels by setting detected artifact intervals to 0
    cfg                        = [];
    cfg.artfctdef.reject       = 'value';
    cfg.artfctdef.eog.artifact = artifact_intervals;
    cfg.artfctdef.value        = 0;
    data_artifact              = ft_rejectartifact(cfg,data_artifact);

end

% apply fieldtrip artifact detection function
%--------------------------------------------------------------------------
for n = 1:num_artifact
    artifact_trl{n}   = detect_artifact_trials(data_artifact,data_artifact.label{n});
end

% save data
%--------------------------------------------------------------------------
artifact_trials.artifact_trials = artifact_trl;
artifact_trials.label           = data_artifact.label';

if save_data
    filename_new = [subjectdata.subjectname '_' files2preproc{m} '_artifact_trials_4ica.mat'];
    save([path2save filesep filename_new],'artifact_trials','-v7.3')
end

% show found artifacts
%--------------------------------------------------------------------------
if check_artifacts
    idx = contains(data_artifact.label,'eog','IgnoreCase',true);
    cfg                        = [];
    cfg.channel                = {'eog'}; % components to be plotted
    cfg.viewmode               = 'vertical';
    cfg.blocksize              = 750;
    cfg.plotevents             = 'no';
    cfg.artfctdef.eog.artifact = vertcat(artifact_trl{idx});
    if additional_suppression
        % if suppression is online show detected artifacts in ecg channel
        cfg.artfctdef.ecg.artifact = artifact_intervals;
    end
    ft_databrowser(cfg,data_artifact);
end

% aks user for input
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

%% functions
%--------------------------------------------------------------------------

function [artifact_intervals] = detect_artifact_trials(data_artifact,channel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performs artifact trial detection based on beforehand specified 
% channeltype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(contains(channel,'ecg','IgnoreCase',true)) % checks if at least one element is nonzero
    cfg                       = [];
    cfg.continuous            = 'yes';
    cfg.artfctdef.ecg.pretim  = 0.25;
    cfg.artfctdef.ecg.psttim  = 0.50-1/1200;
    cfg.artfctdef.ecg.inspect = channel; % channels shown in QRS-lockes average
    cfg.artfctdef.ecg.channel = channel;
    [~, artifact_intervals]   = ft_artifact_ecg(cfg, data_artifact);

elseif any(contains(channel,'eog','IgnoreCase',true)) % checks if at least one element is nonzero
    cfg                       = [];
    cfg.continuous            = 'yes';
    cfg.artfctdef.eog.channel = channel;
    [~, artifact_intervals]   = ft_artifact_eog(cfg, data_artifact);

else
    disp('channel not found')
end  
    
end

%--------------------------------------------------------------------------

function [trl] = detect_threshold_artifacts(cfg,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cfg.threshold in uV
% cfg.pretim    in sec
% cfg.psttim    in sec
% cfg.channel 

% trl: trial defination for artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold = cfg.threshold*10^-6;     % in V
pretim    = data.fsample*cfg.pretim; % in samples
psttim    = data.fsample*cfg.psttim; % in samples

T         = length(data.trial);
trl       = cell(1,T);
offset    = 0;

for t=1:T
    cidx  = strcmp(data.label,cfg.channel);
    trial = data.trial{t}(cidx,:);
    idx   = find(abs(trial)>threshold);
    L     = size(trial,2);
    
    if ~isempty(idx)
        cluster = give_cluster(idx);

        trl{t}    = [cluster(:,1)-pretim,cluster(:,2)+psttim];

        % correct if you are out of the interval
        lidx         = find(trl{t}<1);
        trl{t}(lidx) = length(1);
        hidx         = find(trl{t}>L);
        trl{t}(hidx) = L;

        % correct for offset - see trials as continuous file
        
        trl{t} = trl{t} + offset;
    end
    
    offset = offset + L;
    
end
    
    trl = vertcat(trl{:});       
            
end

%--------------------------------------------------------------------------

function [cluster] = give_cluster(idx)
cluster = [idx(1),0]; 
n       = 1;
clustnr = 1;
    while n<length(idx)
        if idx(n+1)~=idx(n)+1
            cluster(clustnr,2) = idx(n);
            clustnr            = clustnr+1;
            cluster            = [cluster;[idx(n+1),0]];
        end

        if n==(length(idx)-1)
            cluster(end,2) = idx(end);      
        end

        n = n+1;
    end

end
            
       

     
