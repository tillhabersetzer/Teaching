close all
clear 
clc 

%% Import main settings 
%--------------------------------------------------------------------------
addpath('..')
eval('main_settings')

%% Script settings
%--------------------------------------------------------------------------
subidx  = 1;
subject = ['sub-',num2str(subidx,'%02d')];

% artifacts to plot in databrowser
artifacts = {'ECG003','EOG001','EOG002'};

%% Load data
%--------------------------------------------------------------------------
dir2load         = fullfile(settings.path2project,'derivatives',subject,'sensorlevel');
data             = importdata(fullfile(dir2load,[subject,'_sensor_comp_ica.mat']));
data_sensor_comp = data.data_sensor_comp;
coi_artifact     = data.coi_artifact;
clear data

%% 1.) compute correlations
num_artifact = length(coi_artifact);
pearson_corr = cell(1,num_artifact);

for a = 1: num_artifact 
    % number of trials
    num_trials = length(data_sensor_comp{a}.trial);
    % number of components
    num_comp   = length(data_sensor_comp{a}.label)-1;
    % initilaize matrix for correlations
    correlations = zeros(num_comp,num_trials);
    for i = 1:num_trials
        % data for correlattion
        data2corr         = data_sensor_comp{a}.trial{i};
        correlations(:,i) = corr(data2corr(1,:)',data2corr(2:end,:)','Type','Pearson')'; % columswise comparison
    end
    pearson_corr{a} = correlations;
    clear num_trials correlations
end      

%% 2.) plot correlations

% channels of interest for artefact detection
artifact_type = {'ECG','EOG001','EOG002'};
% choose here - values correspond to index in coi_artifact/coi_sensor
%--------------------------------------------------------------------------
artifact = [1,2,3]; % ECG etc.
%--------------------------------------------------------------------------
% find index of elements in arti_comp_data to loop over
ind_a  = contains(coi_artifact,artifact_type(artifact),'IgnoreCase',true)';
a_loop = find(ind_a);

for a = a_loop
    fig_name = [subject ' | Sensor: meg | Artifact: ',coi_artifact{a}];
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
    subtitle(fig_name)
end  

%% 3.) average the components timelocked to the artifact-complex
% only valid for ecg artifacts cause the were gathered timelocked and all
% trials have the same length
cfg   = [];
ind_a = contains(coi_artifact,{'ECG'},'IgnoreCase',true);
a     = find(ind_a);
 
timelock = ft_timelockanalysis(cfg,data_sensor_comp{a});
    
% look at the timelocked/averaged components and compare them with the ECG

fig_name = [subject ' | Sensor: meg | Artifact: ',coi_artifact{a}];
figure('name',fig_name);
title('comparison of timelocked components and ecg')
subplot(2,1,1); plot(timelock.time, timelock.avg(1,:))
subplot(2,1,2); plot(timelock.time, timelock.avg(2:end,:))
subtitle(fig_name)
figure('name',fig_name);
title('comparison of timelocked components and ecg')
subplot(2,1,1); plot(timelock.time, timelock.avg(1,:))
subplot(2,1,2); imagesc(timelock.avg(2:end,:));
subtitle(fig_name)

%% fill out identified bad components
bad_components_meg = [2,21,80];

data = concatenate_components_and_artifacts(dir2load,subject);

%% 4.) plot independent components timecourse with artifact channels 

% artifact trials
trials = importdata(fullfile(dir2load,[subject,'_artifact_trials_4ica.mat']));

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
cfg.mychanscale            = 10^5*ones(length(cfg.mychan),1);                        % scaling of components
ft_databrowser(cfg,data.observe_artifacts);
subtitle([subject ' | component timecourse']);

%% 5.) show components on topoplot
% load components
%----------------

% ica data
components = importdata(fullfile(dir2load,[subject,'_ica.mat']));

% Show components with topoplot
%------------------------------

comp_meg      = components;
cfg           = [];
cfg.component = bad_components_meg;   % specify the component(s) that should be plotted
cfg.layout    = 'Neuromag306mag.lay';  % specify the layout file that should be used for plotting
cfg.comment   = 'no';
figure('name','topoplot meg')
ft_topoplotIC(cfg, comp_meg)
subtitle([subject ' | component topoplot']);

%% show additional components

% time course
%------------
data_sensor = importdata(fullfile(dir2load,[subject,'_data_sensor_4ica.mat']));    

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
subtitle([subject ' | add. component timecourse']);

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

%% apply ICA and compare data

condition = 'upchirps';
datapath  = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',condition,'.fif']);

cfg              = [];
cfg.dataset      = datapath;
cfg.channel      = 'meg'; 
cfg.demean       = 'yes';
cfg.detrend      = 'yes';
cfg.continuous   = 'yes';
cfg.coilaccuracy = 0;       
cfg.bpfilter     = 'yes';
cfg.bpfreq       = [1,150];
cfg.dftfilter    = 'yes';       
cfg.dftfreq      = [50 100 150];
cfg.dftreplace   = 'zero';
data_raw         = ft_preprocessing(cfg);  

% remove the bad components and backproject the data
cfg           = [];
cfg.component = [2,21]; % to be removed component(s)
data_ica      = ft_rejectcomponent(cfg, components, data_raw);

cfg           = [];
cfg.channel   = 'megmag'; 
cfg.viewmode  = 'butterfly';
cfg.blocksize = 10;
ft_databrowser(cfg, data_raw);

ft_databrowser(cfg, data_ica);


%% functions
%--------------------------------------------------------------------------

function [data] = concatenate_components_and_artifacts(dir2load,subject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates data for the detection of artifact related
% independent components. It concatenates the independent components of one 
% sensortype (megmag,megplanar) with the artifact channels
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
data_sensor = importdata(fullfile(dir2load,[subject,'_data_sensor_4ica.mat']));        

% continous artifact data (eog,ecg)
data_artifact = importdata(fullfile(dir2load,[subject,'_data_artifact_4ica.mat']));

% ica unmixing matrices
components = importdata(fullfile(dir2load,[subject,'_ica.mat']));

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
end

%--------------------------------------------------------------------------

function [chan_ind,channel_list_without_artifacts] = find_channel_index(label,bad_components,artifacts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds corresponding indices (chan_ind) in label for 
% specified channel names in bad_components and artifacts. 
% - channel_list_without_artifacts contains the name of the indipendent
%   component channels
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
 % add artifact labels
channel_list                   = [channel_list,artifacts];
% get channel indices
chan_ind                       = find(contains(label,channel_list,'IgnoreCase',true));
end
