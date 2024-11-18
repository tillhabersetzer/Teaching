close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check identified independent components visually
%-------------------------------------------------
% - visualizing identified artifact related independent components
%   together with artifact channels
% - or just independent components topoplot
% - the expected bad components must be entered beforehand in 
%   'Z:\analysis\subject_files\subject_ica_assumptions
% - the detected components must be inserted into the subjectfiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subjects
subjectlist = {'subject24'};

% choose files
% files2preproc = {'stories_maxfilter','olsa_maxfilter'};
files2preproc = {'all_maxfilter'};

% option to generate data for comparison and saving it
generate_data = 1;

% option to save 
save_data = 1;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);

% artifacts to plot in databrowser
artifacts = {'ECG','EOG001','EOG002'};

%% initialize
    
% loop over subjects 
%-------------------
eval('subject_ica_assumptions')
N_subj  = length(subjectlist);
N_files = length(files2preproc);

for s = 1:N_subj
    
% subject selected out of codelist 
subject = subjectlist{s};
eval(subject)

s_fieldnames = fieldnames(icalist);
sidx         = strcmp(s_fieldnames,subjectlist{s});
c_fields     = icalist.(s_fieldnames{sidx});

% loop over selected files
%-------------------------

for m = 1:N_files
    
m_filenames        = fieldnames(c_fields);
midx               = strcmp(m_filenames,files2preproc{m});
bad_components_meg = c_fields.(m_filenames{midx});

if contains(files2preproc{m},'stories')
    path2save = subjectdata.ica_stories;
elseif contains(files2preproc{m},'olsa')
    path2save = subjectdata.ica_olsa;
elseif contains(files2preproc{m},'all')
    path2save = subjectdata.ica_all;
else
    error('Specified path2save does not exist!')
end    

%% 1.) generate structure containing artifacts and independent components if necessary
if generate_data
    data = concatenate_components_and_artifacts(subjectdata,files2preproc{m},path2save);
    if save_data
        filename_new = [subjectdata.subjectname '_' files2preproc{m} '_observe_artifacts_4ica.mat'];
        save([path2save filesep filename_new],'data','-v7.3')
    end 
end

if ~generate_data
    filename_new = [subjectdata.subjectname '_' files2preproc{m} '_observe_artifacts_4ica.mat'];
    if ~isfile([path2save filesep filename_new])
        error('Please generate observe_artifacts_4ica.')
    end
    data = importdata([path2save filesep filename_new]);
end 

%% 2.) plot independent components timecourse with artifact channels 

% artifact trials
filename_trials = [subjectdata.subjectname '_' files2preproc{m} '_artifact_trials_4ica.mat'];
trials          = importdata([path2save filesep filename_trials]);

% find correct index positions for artifact trials
ind_eog = contains(trials.label,'eog','IgnoreCase',true);
ind_ecg = contains(trials.label,'ecg','IgnoreCase',true);
trl_eog = trials.artifact_trials(ind_eog);
trl_eog = [trl_eog{1};trl_eog{2}]; % both eog channels together
trl_ecg = trials.artifact_trials{ind_ecg};

% meg
%--------------------------------------------------------------------------
label = data.observe_artifacts.label;

cfg                        = [];
cfg.continuous             = 'yes';
cfg.viewmode               = 'vertical';
cfg.blocksize              = 50;        % length of data to display, in sec
% cfg.artfctdef.ecg.artifact = trl_ecg;  % easy to see
cfg.artfctdef.eog.artifact = trl_eog;   % here you can mark artifacts
cfg.eogscale               = 1;         % eog scaling es reference
cfg.ecgscale               = 0.3;
[cfg.channel,~]            = find_channel_index(label,bad_components_meg,artifacts); % index of all channels
% treat components separate to scale them appropriate
[~,cfg.mychan]             = find_channel_index(label,bad_components_meg,{});        % names of component channels 
cfg.mychanscale            = 10^5*ones(length(cfg.mychan),1);                          % scaling of components
ft_databrowser(cfg,data.observe_artifacts);
suptitle([subjectlist{s} ' | ' files2preproc{m} ' | component timecourse']);

%% show components on topoplot
% load components
%----------------

% ica data
filename_ica = [subjectdata.subjectname '_' files2preproc{m} '_ica_comp.mat'];
components   = importdata([path2save filesep filename_ica]);

% Show components with topoplot
%------------------------------

comp_meg      = components;
cfg           = [];
cfg.component = bad_components_meg;   % specify the component(s) that should be plotted
cfg.layout    = 'Neuromag306mag.lay';  % specify the layout file that should be used for plotting
cfg.comment   = 'no';
figure('name','topoplot meg')
ft_topoplotIC(cfg, comp_meg)
suptitle([subjectlist{s} ' | ' files2preproc{m} ' | component topoplot']);


%% show additional components

% time course
%------------
filename_data = [subjectdata.subjectname '_' files2preproc{m} '_data_sensor_4ica.mat'];
data_sensor   = importdata([path2save filesep filename_data]);  

cfg           = [];
cfg.channel   = {'meg'};
cfg.unmixing  = components.unmixing;
cfg.topolabel = components.topolabel;
comp          = ft_componentanalysis(cfg,data_sensor);

cfg            = [];
cfg.continuous = 'yes';
cfg.blocksize  = 100;   
cfg.layout     = 'Neuromag306mag.lay'; % specify the layout file that should be used for plotting
cfg.viewmode   = 'component';
ft_databrowser(cfg,comp)
suptitle([subjectlist{s} ' | ' files2preproc{m} ' | add. component timecourse']);

% topoplot
%---------
% define additional components, just for sanity
extra_components = 1:length(comp_meg.label);

comp_meg      = components;
cfg           = [];
cfg.component = extra_components;   % specify the component(s) that should be plotted
cfg.layout    = 'Neuromag306mag.lay';  % specify the layout file that should be used for plotting
cfg.comment   = 'no';
figure('name','topoplot meg')
ft_topoplotIC(cfg, comp_meg)
suptitle([subjectlist{s} ' | ' files2preproc{m} ' | add. component topoplot']);

% aks user for input
response = input('Next file? [y]/[n]:','s');
if strcmpi(response,'n')
    break
end

close all

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

function [data] = concatenate_components_and_artifacts(subjectdata,filename,path2save)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate data for artifact detection
% concatenates the independent components of one sensortype 
% (megmag,megplanar) with the artifact channels
% this enables the option to compare the independent components with the
% artifact channels
%
% - data_sensor: contains data from all sensor channels
%   (megmag,megplanar)
% - data_artifact: contains data of artifact channels
% - components: contains topolabels and unmixing matrices for reach
%   sensortype
% - save_data = 1 or 0, enables option to save generated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
%--------------------------------------------------------------------------
% continuous data (meg)
filename_data = [subjectdata.subjectname '_' filename '_data_sensor_4ica.mat'];
data_sensor   = importdata([path2save filesep filename_data]);        

% continous artifact data (eog,ecg)
filename_artifact = [subjectdata.subjectname '_' filename '_data_artifact_4ica.mat'];
data_artifact     = importdata([path2save filesep filename_artifact]);

% ica unmixing matrices
filename_ica = [subjectdata.subjectname '_' filename '_ica_comp.mat'];
components   = importdata([path2save filesep filename_ica]);

% concatenate data  
% decompose data into independent components
cfg           = [];
cfg.channel   = {'meg'};
cfg.unmixing  = components.unmixing;
cfg.topolabel = components.topolabel;
comp_data     = ft_componentanalysis(cfg,data_sensor);

% append artifact data and independet components
cfg               = [];
observe_artifacts = ft_appenddata(cfg,data_artifact,comp_data); 
  
data.observe_artifacts = observe_artifacts;
data.label             = {'meg'};
data.triallabel        = data_sensor.triallabel;
end

%--------------------------------------------------------------------------

function [chan_ind,channel_list_without_artifacts] = find_channel_index(label,bad_components,artifacts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finds index (chan_ind) in label for specified channels in bad_components
% and artifacts. channel_list_without_artifacts contains the complete name
% of the desired independent components to look at
% - artifacts must be a cell array like {'eog','ecg'}
% - bad_components is an array containing detected bad components like
%   [3,10] for component 3 and component 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name1        = 'component';
channel_list = {};

for n = bad_components
    if n < 10
        name2 = ['00' num2str(n)];
    elseif (9 < n) && (n < 100)
        name2 = ['0' num2str(n)];
    else
        name2 = num2str(n);
    end
   channel_list = [channel_list,{[name1,name2]}];
end

channel_list_without_artifacts = channel_list;
channel_list = [channel_list,artifacts];
check_arti   = contains(label,channel_list,'IgnoreCase',true);
chan_ind     = find(check_arti);

end





         