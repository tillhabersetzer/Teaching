close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
for i = 3:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjects{i-2} = [subject,num2str(i)]; 
end

check = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'));
% addpath for preprocessing function
addpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));


% loop over subjects 
%-------------------
N_subj       = length(subjects);
sources      = cell(1,N_subj);
sources_pre  = cell(1,N_subj);
sources_pst  = cell(1,N_subj);
sources_diff = cell(1,N_subj);
parcel       = cell(1,N_subj);
info         = cell(1,N_subj);

for s = 1:N_subj
    eval(subjects{s})
    source              = importdata(fullfile(subjectdata.chirps_mne,[subjectdata.subjectname,'_mne_sources.mat']));
    sources{s}          = source.source;
    
    sources_pre{s}      = source.source_pre;
    sources_pst{s}      = source.source_pst;
    
    sources{s}.pos      = sources{1}.pos;
    sources_pre{s}.pos  = sources{1}.pos;
    sources_pst{s}.pos  = sources{1}.pos;
    % change time axis
%     sources_pre{s}.time = sources_pst{1}.time;
    sources_pst{s}.time = sources_pre{1}.time;
    %sources_pre{s} = rmfield(sources_pre{s},'inside');
    %sources_pst{s} = rmfield(sources_pst{s},'inside');
    %sources_pre{s}.avg = rmfield(sources_pre{s}.avg,'label');
    %sources_pst{s}.avg = rmfield(sources_pst{s}.avg,'label');
    
    cfg            = [];
    % project the source to its strongest orientation, i.e. the direction that explains most of the source variance. 
    % That is equivalent to taking the largest eigenvector of the source
    % timeseries.
    cfg.projectmom = 'yes';
    sources    {s} = ft_sourcedescriptives(cfg,sources{s});
    sources_pre{s} = ft_sourcedescriptives(cfg,sources_pre{s});
    sources_pst{s} = ft_sourcedescriptives(cfg,sources_pst{s});
    
    cfg                 = [];
    cfg.parameter       = 'pow';
    %cfg.operation       = '((x1-x2)./max(abs(x1-x2)))'; %percentage change from baseline
    cfg.operation       = 'subtract';
    sources_diff{s}     = ft_math(cfg,sources_pst{s},sources_pre{s});
    sources_diff{s}.pow = sources_diff{s}.pow./max(abs(sources_diff{s}.pow),[],'all');
    info{s}             = source.info;   
end

% calculate grand average  
%------------------------
cfg           = [];
%cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'pow';
grandavg      = ft_sourcegrandaverage(cfg, sources{:});
grandavg_pre  = ft_sourcegrandaverage(cfg, sources_pre{:});
grandavg_pst  = ft_sourcegrandaverage(cfg, sources_pst{:});
grandavg_diff = ft_sourcegrandaverage(cfg, sources_diff{:});

grandavg.tri      = sources{1,1}.tri;
grandavg_pre.tri  = sources{1,1}.tri;
grandavg_pst.tri  = sources{1,1}.tri;
grandavg_diff.tri = sources{1,1}.tri;


% plot data
%----------
cfg              = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,sources_diff{22});

cfg              = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,grandavg);

cfg              = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,grandavg_pre);

cfg              = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,grandavg_pst);

cfg              = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,grandavg_diff);

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'pow';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.latency      = [0.05,0.15];
cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,grandavg_diff);

cfg            = [];
cfg.projectmom = 'yes';
sd             = ft_sourcedescriptives(cfg,sources{1});

%% statistics - DOES NOT WORK - SEEMS TO HAVE SEVERE PROBLEMS WITH SURFACE DATA
%-----------
cfg                         = [];
cfg.parameter               = 'pow';
% cfg.dim                     = sources_pre{1}.dim;
cfg.method                  = 'montecarlo';
cfg.statistic               = 'ft_statfun_depsamplesT';
cfg.correctm                = 'cluster';
cfg.clusteralpha            = 0.05;
%cfg.clusterstatistic        = 'maxsum';
cfg.clusterstatistic        = 'max';
cfg.tail                    = 0;
cfg.clustertail             = 0;
cfg.alpha                   = 0.025;
cfg.numrandomization        = 10;
design                      = zeros(2,2*N_subj);
design(1,1:N_subj)          = 1;
design(1,N_subj+1:2*N_subj) = 2;
design(2,1:N_subj)          = [1:N_subj];
design(2,N_subj+1:2*N_subj) = [1:N_subj];
cfg.design                  = design;
cfg.ivar                    = 1; % row of design matrix that contains independent variable 
cfg.uvar                    = 2; % row of design matrix that contains unit of observation
cfg.tri                     = sources{1,1}.tri;
cfg.inside                  = sources{1,1}.inside;

stat                        = ft_sourcestatistics(cfg,sources_pst{:},sources_pre{:});


% run statistics over subjects %
cfg                  = [];
%cfg.channel          = 'all';
cfg.neighbours       = []; % no channel neighbours, only time
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.parameter        = 'pow';
cfg.correctm         = 'cluster';
cfg.correcttail      = 'prob';
cfg.numrandomization = 10;
cfg.alpha            = 0.05; % note that this only implies single-sided testing

cfg.design(1,:) = [1:N_subj 1:N_subj];
cfg.design(2,:) = [ones(1,N_subj)*1 ones(1,N_subj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat = ft_sourcestatistics(cfg, sources_pst{:}, sources_pre{:});





    
%% Clean up
rmpath(fullfile('Z:','analysis','subject_files'));
rmpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));