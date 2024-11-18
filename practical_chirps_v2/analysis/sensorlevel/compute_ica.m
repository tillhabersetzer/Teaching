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

conditions = settings.conditions;
C          = length(conditions);


% these cells are saved for further processing
data_sensor   = cell(1,C);
data_artifact = cell(1,C);
% needed for ica on all files
data_epoched  = cell(1,C);

%% Computations
%--------------------------------------------------------------------------
for cidx=1:C

    datapath = fullfile(settings.path2project,'rawdata',subject,'meg',[subject,'_task-',conditions{cidx},'.fif']);

    %% 1.) load and filter continuous data
    cfg              = [];
    cfg.dataset      = datapath;
    % cfg.channel      = [{'meg','eog','ecg'},strcat('-',subjectdata.badchannels)]; 
    cfg.channel      = {'meg','eog','ecg'}; 
    cfg.demean       = 'yes';
    cfg.detrend      = 'yes';
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0;       
    cfg.bpfilter     = 'yes';
    cfg.bpfreq       = [1,150];
    cfg.dftfilter    = 'yes';       
    cfg.dftfreq      = [50 100 150];
    cfg.dftreplace   = 'zero';
    data_filtered    = ft_preprocessing(cfg);  

    %% 2.) separate the data
    cfg                 = [];
    cfg.channel         = {'ecg','eog'};
    data_artifact{cidx} = ft_selectdata(cfg, data_filtered); 
    cfg.channel         = 'meg';
    data_sensor{cidx}   = ft_selectdata(cfg, data_filtered); 

    %% 3.) define segment for ica
    cfg                      = [];
    cfg.dataset              = datapath;
    cfg.trialfun             = 'ft_trialfun_general';
    cfg.trialdef.triallength = 5;
    cfg                      = ft_definetrial(cfg);
    trl                      = cfg.trl;

    %% 4.) epoch data
    cfg                = [];
    % multiple segmentes between triggers
    cfg.trl            = trl;
    data_epoched{cidx} = ft_redefinetrial(cfg,data_sensor{cidx});
    
    cfg                 = [];
    % entire segments between triggers
    cfg.trl             = [trl(1,1),trl(end,2),0];
    data_artifact{cidx} = ft_redefinetrial(cfg,data_artifact{cidx});
    
    cfg               = [];
    % entire segments between triggers
    cfg.trl           = [trl(1,1),trl(end,2),0];
    data_sensor{cidx} = ft_redefinetrial(cfg,data_sensor{cidx});

    %% 5.) downsample data to reduce memory load
    cfg                 = [];
    cfg.resamplefs      = 400;
    cfg.detrend         = 'no';
    data_epoched{cidx}  = ft_resampledata(cfg,data_epoched{cidx}); 
    data_artifact{cidx} = ft_resampledata(cfg,data_artifact{cidx});
    data_sensor{cidx}   = ft_resampledata(cfg,data_sensor{cidx});
 
    %% 6.) clean ica data - reject heavily contaminated epochs    

    cfg                        = [];
    cfg.artfctdef.jump.channel = 'meg'; 
    [~, artifact_jump]         = ft_artifact_jump(cfg, data_epoched{cidx});
    
    cfg                        = [];
    cfg.artfctdef.clip.channel = 'all'; 
    [~, artifact_clip]         = ft_artifact_clip(cfg, data_epoched{cidx});
    
    cfg                              = [];
    cfg.artfctdef.zvalue.interactive = 'no';
    cfg.artfctdef.zvalue.channel     = 'all'; 
    cfg.artfctdef.zvalue.cutoff      = 15;
    [~, artifact_zvalue]             = ft_artifact_zvalue(cfg, data_epoched{cidx});
    
    cfg                           = [];
    % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
    cfg.artfctdef.reject          = 'partial';
    cfg.artfctdef.jump.artifact   = artifact_jump;
    cfg.artfctdef.clip.artifact   = artifact_clip;
    cfg.artfctdef.zvalue.artifact = artifact_zvalue;
    data_epoched{cidx}            = ft_rejectartifact(cfg,data_epoched{cidx});   

end

%% 7.) append data
hdr_sensor   = data_sensor{1}.hdr;
hdr_artifact = data_artifact{1}.hdr;

cfg                = [];
cfg.keepsampleinfo = 'no';
data_sensor        = ft_appenddata(cfg,data_sensor{:});
data_artifact      = ft_appenddata(cfg,data_artifact{:});

% add additional information to data
data_sensor.hdr   = hdr_sensor;
data_artifact.hdr = hdr_artifact;

% concatenate data of all trials for ica
cfg                = [];
cfg.keepsampleinfo = 'no';
data_epoched       = ft_appenddata(cfg,data_epoched{:});

% Save data
%----------
% make folder for data
dir2save = fullfile(settings.path2project,'derivatives',subject,'sensorlevel');
if ~exist(dir2save,'dir')
    mkdir(dir2save)
end
save(fullfile(dir2save,[subject,'_data_artifact_4ica.mat']),'data_artifact','-v7.3')
save(fullfile(dir2save,[subject,'_data_sensor_4ica.mat']),'data_sensor','-v7.3')   

%% 7.) ica on different sensortypes
% perform ica
cfg                 = [];
% perform pca for maxfiltered meg data
cfg.runica.pca      = 100; 
% cfg.runica.pca      = r; % rank of data
cfg.runica.sphering = 'on'; % is default (for different sensortypes, variances are equalized)
cfg.method          = 'runica';
cfg.channel         = 'meg';
components          = ft_componentanalysis(cfg,data_epoched);  

% show components
cfg           = [];
cfg.component = [1:30];       % specify the component(s) that should be plotted
cfg.layout    = 'neuromag306mag.lay'; % specify the layout file that should be used for plotting
ft_topoplotIC(cfg, components)

cfg          = [];
cfg.channel  = [1:30]; % components to be plotted
cfg.viewmode = 'component';
cfg.layout   = 'neuromag306mag.lay'; % specify the layout file that should be used for plotting
ft_databrowser(cfg, components)

save(fullfile(dir2save,[subject,'_ica.mat']),'components')

%% 8.) extract artifact epochs
num_artifact = length(data_artifact.label);
artifact_trl = cell(1,num_artifact);

for n = 1:num_artifact
    artifact_trl{n} = detect_artifact_trials(data_artifact,data_artifact.label{n});
end

artifact_trials.artifact_trials = artifact_trl;
artifact_trials.label           = data_artifact.label';
save(fullfile(dir2save,[subject,'_artifact_trials_4ica.mat']),'artifact_trials','-v7.3')
clear artifact_trials

%% 9.) Separate channel by beforehand calculated artifact trials
coi_artifact = data_artifact.label;
% artifact channels
data_artifact_trl = cell(1,num_artifact);
% sensor channels
data_sensor_trl   = cell(1,num_artifact);

for a = 1: num_artifact  
    % only channel of interest
    cfg                  = [];
    cfg.channel          = coi_artifact{a}; 
    data_artifact_trl{a} = ft_selectdata(cfg, data_artifact); 
    
    artifact_trials      = artifact_trl{a};
    cfg                  = [];
    cfg.trl              = [artifact_trials, zeros(size(artifact_trials,1),1)];
    data_artifact_trl{a} = ft_redefinetrial(cfg,data_artifact_trl{a});
         
    % only data segment around artifact trials
    cfg                = [];
    cfg.trl            = [artifact_trials, zeros(size(artifact_trials,1),1)];
    data_sensor_trl{a} = ft_redefinetrial(cfg,data_sensor);
end

%% 10.) decompose the artifact-locked datasegments into components, using the previously found (un)mixing matrix
data_sensor_comp = cell(1,num_artifact);

cfg           = [];
cfg.unmixing  = components.unmixing;
cfg.topolabel = components.topolabel;
for a = 1: num_artifact
  data_sensor_comp{a} = ft_componentanalysis(cfg, data_sensor_trl{a});
end  

%% 11.) append artifact-channels to the sensor data
cfg                = [];
cfg.keepsampleinfo = 'no'; % because data originates from different datafiles

for a = 1:num_artifact
    data_sensor_comp{a} = ft_appenddata(cfg, data_artifact_trl{a},data_sensor_comp{a});
end

% save data
save(fullfile(dir2save,[subject,'_sensor_comp_ica.mat']),'data_sensor_comp','coi_artifact','-v7.3'); 

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
