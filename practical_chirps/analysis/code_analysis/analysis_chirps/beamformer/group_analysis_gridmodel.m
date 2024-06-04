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

% load template information
%--------------------------
template_mri   = ft_read_mri(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ... % mm
                            'template\anatomy\single_subj_T1.nii']);
template_mri_2 = ft_read_mri('Z:/Software/MEG_EEG_Toolboxen/fieldtrip-20191127/external/spm8/templates/T1.nii');
template_grid  = importdata(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ...
                            'template\sourcemodel\standard_sourcemodel3d10mm.mat']);
template_grid  = ft_convert_units(template_grid,'mm'); 
atlas          = ft_read_atlas('Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\template\atlas\aal\ROI_MNI_V4.nii'); % mm

% loop over subjects 
%-------------------
N_subj      = length(subjects);
sources     = cell(1,N_subj);
sources_pre = cell(1,N_subj);
sources_pst = cell(1,N_subj);
info        = cell(1,N_subj);
kappas      = zeros(N_subj,4);

for s = 1:N_subj
    eval(subjects{s})
    source          = importdata(fullfile(subjectdata.chirps_beamformer,...
                                     [subjectdata.subjectname,'_beamformer_gridmodel.mat']));
    cfg             = [];
    cfg.operation   = 'abs';
    cfg.parameter   = 'mom';
    sourcepreabsmom = ft_math(cfg,source.sourcepre);
    sourcepstabsmom = ft_math(cfg,source.sourcepst);
    
    cfg             = [];
    cfg.parameter   = 'mom';
    %cfg.operation   = '((x1-x2)./x2)*100'; %percentage change from baseline
    cfg.operation   = 'subtract';
    sources{s}      = ft_math(cfg,sourcepstabsmom,sourcepreabsmom);
    
    sources{s}.pos  = template_grid.pos;
    info{s}         = source.info;  
    kappas(s,1:2)   = info{s}.kappa_noise;
    kappas(s,3:4)   = info{s}.kappa_data;
    clear source
end
kappas = array2table(kappas,'VariableNames',{'noise-mag','noise-grad','data-mag','data-grad'});

% calculate grand average  
%------------------------
cfg           = [];
cfg.latency   = 'all';
cfg.parameter = 'mom';
grandavg      = ft_sourcegrandaverage(cfg,sources{:});

% generate matrix for moments!!!
I       = length(grandavg.inside);
moments = NaN(I,length(grandavg.time));
for i = find(grandavg.inside)'
        moments(i,:) = grandavg.mom{i};   
end
grandavg.mom      = moments;

%% average without parcellation
%------------------------------

% interpolate onto mri
%---------------------
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'mom';
cfg.interpmethod = 'nearest';
grandavg_int     = ft_sourceinterpolate(cfg,grandavg,template_mri);
grandavg_int.coordsys = 'mni';

% Plot results
%-------------
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.atlas        = atlas;
%cfg.latency      = [0.05,0.15];
%cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,grandavg_int);

% project on brain surface
%-------------------------
grandavg_int_surf     = grandavg_int;
idx                   = dsearchn(grandavg_int_surf.time',[0.05,0.25]'); % timewindow for average
grandavg_int_surf.mom = mean(grandavg_int_surf.mom(:,idx(1):idx(2)),2);

cfg = [];
cfg.method       = 'surface';
cfg.atlas        = atlas;
cfg.funparameter = 'mom';
%cfg.funcolorlim  = [-30 30];
cfg.funcolormap  = 'jet';
cfg.projmethod   = 'nearest';
cfg.surfinflated = 'surface_inflated_both_caret.mat';
cfg.projthresh   = 0.6; % in percent - adapt threshold!!!
cfg.camlight     = 'no';
ft_sourceplot(cfg, grandavg_int_surf);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull    
    
%% Average with parcellation
%---------------------------
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'mom';
cfg.interpmethod = 'nearest';
grandavg_int_2   = ft_sourceinterpolate(cfg, grandavg, template_mri_2);

cfg    = [];
parcel = ft_sourceparcellate(cfg,grandavg_int_2, atlas);

% create a dummy structure where we identify the moment values per voxel
% but remember - moments is a timeseries
idx = dsearchn(grandavg_int_2.time',[0.05,0.15]'); % timewindow for average
dummy = atlas;
for i=1:size(parcel.mom,1)
      dummy.tissue(find(dummy.tissue==i)) = mean(parcel.mom(i,idx(1):idx(2)),2); % arbitrary time value or mean??????
end

% Plot results
%-------------
grandavg_int_2.parcel   = dummy.tissue;
grandavg_int_2.coordsys = 'mni';
cfg                     = [];
cfg.method              = 'ortho';
cfg.funparameter        = 'parcel';
cfg.funcolormap         = 'jet';
cfg.renderer            = 'zbuffer';
cfg.location            = [-42 -20 6];
cfg.atlas               = atlas; % works with atlas!!!
%cfg.funcolorlim        = [-30 30];
ft_sourceplot(cfg,grandavg_int_2);

% project on brain surface
%-------------------------
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
%cfg.funcolorlim    = [-30 30];
cfg.funcolormap    = 'jet';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projthresh     = 0.6; % adapt threshold!!!
cfg.camlight       = 'no';
ft_sourceplot(cfg,grandavg_int_2);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull    

%% statistical thresholding
% load data again
%----------------
for s = 1:N_subj
    eval(subjects{s})
    source         = importdata(fullfile(subjectdata.chirps_beamformer,...
                     [subjectdata.subjectname,'_beamformer_gridmodel.mat']));
    sources_pre{s} = source.sourcepre;
    sources_pst{s} = source.sourcepst;
    % generate matrix for moments!!!
    I           = length(sources_pre{s}.inside);
%     moments_pre = NaN(I,length(sources_pre{s}.time));
%     moments_pst = NaN(I,length(sources_pre{s}.time));
    moments_pre = NaN(I,1);
    moments_pst = NaN(I,1);
    idx         = dsearchn(sources_pre{s}.time',[0.05,0.15]'); % timewindow for average

    for i = find(sources_pre{s}.inside)'
%             moments_pre(i,:) = sources_pre{s}.avg.mom{i}; 
%             moments_pst(i,:) = sources_pst{s}.avg.mom{i};
            moments_pre(i,:) = mean(sources_pre{s}.avg.mom{i}(idx(1):idx(2)).^2); 
            moments_pst(i,:) = mean(sources_pst{s}.avg.mom{i}(idx(1):idx(2)).^2);
    end
    sources_pre{s}.avg.mom = moments_pre;
    sources_pst{s}.avg.mom = moments_pst;
    sources_pre{s}.pos     = template_grid.pos;
    sources_pst{s}.pos     = template_grid.pos;   
end 

% calculate grand average  
%------------------------
cfg           = [];
cfg.latency   = 'all';
cfg.parameter = 'mom';
grandavg_pre  = ft_sourcegrandaverage(cfg,sources_pre{:});
grandavg_pst  = ft_sourcegrandaverage(cfg,sources_pst{:});

% interpolate onto mri
%---------------------
cfg                       = [];
cfg.voxelcoord            = 'no';
cfg.parameter             = 'mom';
cfg.interpmethod          = 'nearest';
grandavg_pre_int          = ft_sourceinterpolate(cfg,grandavg_pre,template_mri);
grandavg_pre_int.coordsys = 'mni';
grandavg_pst_int          = ft_sourceinterpolate(cfg,grandavg_pst,template_mri);
grandavg_pst_int.coordsys = 'mni';

% Plot results
%-------------
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.atlas        = atlas;
%cfg.latency      = [0.05,0.15];
%cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,grandavg_pre_int);
ft_sourceplot(cfg,grandavg_pst_int);

% statistics
%-----------
cfg                         = [];
cfg.parameter               = 'mom';
cfg.dim                     = sources_pre{1}.dim;
cfg.method                  = 'montecarlo';
cfg.statistic               = 'ft_statfun_depsamplesT';
cfg.correctm                = 'cluster';
cfg.clusteralpha            = 0.05;
%cfg.clusterstatistic        = 'maxsum';
cfg.clusterstatistic        = 'max';
cfg.tail                    = 0;
cfg.clustertail             = 0;
cfg.alpha                   = 0.025;
cfg.numrandomization        = 1000;
design                      = zeros(2,2*N_subj);
design(1,1:N_subj)          = 1;
design(1,N_subj+1:2*N_subj) = 2;
design(2,1:N_subj)          = [1:N_subj];
design(2,N_subj+1:2*N_subj) = [1:N_subj];
cfg.design                  = design;
cfg.ivar                    = 1; % row of design matrix that contains independent variable 
cfg.uvar                    = 2; % row of design matrix that contains unit of observation
stat                        = ft_sourcestatistics(cfg,sources_pst{:},sources_pre{:});

% without parcellation
%--------------------------------------------------------------------------

% interpolation of significant voxels
%------------------------------------
%stat.inside       = template_grid.inside;
cfg               = [];
cfg.voxelcoord    = 'no';
cfg.parameter     = 'stat';
cfg.interpmethod  = 'nearest';
stat_int          = ft_sourceinterpolate(cfg, stat, template_mri); % interpolate result
cfg.parameter     = 'mask';
mask_int          = ft_sourceinterpolate(cfg, stat, template_mri); % interpolate mask
stat_int.mask     = mask_int.mask;
stat_int.coordsys = 'mni';

% plot result
%------------
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.maskparameter = 'mask';
cfg.atlas         = atlas;
cfg.location      = 'max';
cfg.funcolorlim   = [-5 5];
cfg.funcolormap   = 'jet';
ft_sourceplot(cfg,stat_int);

% project on brain surface - not possible for timeseries
%-------------------------------------------------------
cfg               = [];
cfg.method        = 'surface';
cfg.atlas         = atlas;
cfg.funparameter  = 'stat';
cfg.maskparameter = 'mask';
%cfg.funcolorlim  = [-30 30];
cfg.funcolormap   = 'jet';
cfg.projmethod    = 'nearest';
cfg.surfinflated  = 'surface_inflated_both_caret.mat';
cfg.projthresh    = 0.8; % in percent - adapt threshold!!!
cfg.camlight      = 'no';
ft_sourceplot(cfg, stat_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull    

%% with parcellation - not possible for timeseries
%-------------------------------------------------

% interpolation of significant voxels
%------------------------------------
%stat.inside         = template_grid.inside;
cfg                 = [];
cfg.voxelcoord      = 'no';
cfg.parameter       = 'stat';
cfg.interpmethod    = 'nearest';
stat_int_2          = ft_sourceinterpolate(cfg, stat, template_mri_2); % interpolate result
% cfg.parameter       = 'mask';
% mask_int_2          = ft_sourceinterpolate(cfg, stat, template_mri_2); % interpolate mask
% stat_int_2.mask     = mask_int_2.mask;
stat_int_2.coordsys = 'mni';

cfg        = [];
parcel     = ft_sourceparcellate(cfg, stat_int_2, atlas);
% parcelmask = ft_sourceparcellate(cfg, mask_int_2, atlas);

% create dummy struct
dummy     = atlas;
dummymask = atlas;
for i=1:length(parcel.stat)
      dummy.tissue(find(dummy.tissue==i))        = parcel.stat(i);
%       dummymask.tissue(find(dummymask.tissue==i))= (parcelmask.mask(i));
end

% plot result
%------------
stat_int_2.parcel   = dummy.tissue;
stat_int_2.coordsys = 'mni';
%stat_int.mask     = dummymask.tissue;
cfg               = [];
cfg.method        = 'slice';
cfg.funparameter  = 'parcel';
cfg.funcolormap   = 'jet';
%cfg.maskparameter = 'mask';
cfg.renderer      = 'zbuffer';
cfg.funcolorlim   = [-5 5];
cfg.atlas         = atlas;
ft_sourceplot(cfg,stat_int_2);



    
%% Clean up
rmpath(fullfile('Z:','analysis','subject_files'));
rmpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));