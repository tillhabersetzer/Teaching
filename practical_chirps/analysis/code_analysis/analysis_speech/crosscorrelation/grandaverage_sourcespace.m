close all; clear all; clc;

% Settings
%--------------------------------------------------------------------------
% choose subject 
% subjectlist = {'subject04'};
for i = 2:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjectlist{i-1} = [subject,num2str(i)]; 
end

% choose files
files2preproc = 'stories_maxfilter';

% choose envelope type for crosscorrelation
% envelopetype = 'onset_envelope';
envelopetype = 'envelope';

% load data with performed ica (1) or without (0)
ica_on = 1; 

% load data with whitening (1) or without (0)
% not necessary for correlation in sensorspace! cross correlation is
% normalized
whitening = 1; 

% load data with additional zscoring of all trials
zscoring = 1;

% timewindow
timewin = [0.03,0.12;0.12,0.19];
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'))

if ica_on
    add = '_ica';
else 
    add = '';
end

if zscoring
    add2 = '_zscored';
else 
    add2 = '';
end

new_dir = fullfile('Z:','analysis','generated_data','grouplevel','speech',...
          'crosscorrelation',envelopetype,'sourcelevel');
if ~exist(new_dir, 'dir')
   mkdir(new_dir)
end

% load template information for beamformer
%-----------------------------------------
template_mri   = ft_read_mri(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ... % mm
                            'template\anatomy\single_subj_T1.nii']);
template_mri_2 = ft_read_mri('Z:/Software/MEG_EEG_Toolboxen/fieldtrip-20191127/external/spm8/templates/T1.nii');
template_grid  = importdata(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ...
                            'template\sourcemodel\standard_sourcemodel3d10mm.mat']);
template_grid  = ft_convert_units(template_grid,'mm'); 
atlas          = ft_read_atlas('Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\template\atlas\aal\ROI_MNI_V4.nii'); % mm
template_surf  = importdata(fullfile('Z:','analysis','generated_data','single_subj_T1_1mm','sourcemodel',...
                                     'single_subj_T1_1mm_sourcemodel_4k.mat'));

%% load data
%--------------------------------------------------------------------------
N_subj      = length(subjectlist);
avg         = cell(1,N_subj);
avg_shuffle = cell(1,N_subj);
info        = cell(1,N_subj);

for s = 1:N_subj
    eval(subjectlist{s})
    if ica_on
       filename_new = [subjectdata.subjectname '_' files2preproc '_crosscorr_ica'];
    else
       filename_new = [subjectdata.subjectname '_' files2preproc '_crosscorr_noica'];
    end
    if whitening
        filename_new = [filename_new '_whitened'];
    end
    if zscoring
        filename_new = [filename_new '_zscored'];
    end
    % data
    %-----
    data               = importdata(fullfile(subjectdata.speech,'crosscorrelation',...
                         envelopetype,'sourcelevel',[filename_new,'.mat']));
    avg{s}             = data.source_avg;
    avg_shuffle{s}     = data.source_avg_shuffle;
    avg{s}.pos         = template_grid.pos;
    avg_shuffle{s}.pos = template_grid.pos;
    info{s}            = data.info;  
    clear data
    disp([subjectlist{s},' loaded.'])
end

%% grandaverage
%--------------------------------------------------------------------------
% beamformer
%-----------
cfg              = [];
cfg.latency      = 'all';
cfg.parameter    = 'mom';
grandavg         = ft_sourcegrandaverage(cfg,avg{:});
grandavg_shuffle = ft_sourcegrandaverage(cfg,avg_shuffle{:});

cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'mom';
raweffect     = ft_math(cfg,grandavg,grandavg_shuffle);

%% visualize results
%-------------------
% interpolate onto mri
cfg                            = [];
cfg.voxelcoord                 = 'no';
cfg.parameter                  = 'mom';
cfg.interpmethod               = 'nearest';
% grandavg_int                   = ft_sourceinterpolate(cfg,grandavg,template_mri);
grandavg_shuffle_int           = ft_sourceinterpolate(cfg,grandavg_shuffle,template_mri);
raweffect_int                  = ft_sourceinterpolate(cfg,raweffect,template_mri);
% grandavg_int.coordsys          = 'mni';
grandavg_shuffle_int.coordsys  = 'mni';
raweffect_int.coordsys         = 'mni';

% Plot results
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.atlas        = atlas;
% cfg.latency     = timewin;
% cfg.avgovertime = 'yes';
% ft_sourceplot(cfg,grandavg_int);
ft_sourceplot(cfg,grandavg_shuffle_int);
ft_sourceplot(cfg,raweffect_int);

%--------------------------- with parcellation-----------------------------
cfg                    = [];
cfg.voxelcoord         = 'no';
cfg.parameter          = 'mom';
cfg.interpmethod       = 'nearest';
raweffect_int          = ft_sourceinterpolate(cfg,raweffect,template_mri_2);
raweffect_int.coordsys = 'mni';

cfg    = [];
parcel = ft_sourceparcellate(cfg,raweffect_int,atlas);

% create a dummy structure where we identify the moment values per voxel
% but remember - moments is a timeseries
dummy1 = atlas;
dummy2 = atlas;
idx1   = dsearchn(raweffect_int.time',timewin(1,:)');
idx2   = dsearchn(raweffect_int.time',timewin(2,:)');
for i=1:size(parcel.mom,1)
      dummy1.tissue(find(dummy1.tissue==i)) = mean(parcel.mom(i,idx1(1):idx1(2))); 
      dummy2.tissue(find(dummy2.tissue==i)) = mean(parcel.mom(i,idx2(1):idx2(2))); 
end

cfg                  = [];
cfg.method           = 'ortho';
cfg.funparameter     = 'parcel';
cfg.funcolormap      = 'jet';
cfg.renderer         = 'zbuffer';
cfg.location         = [-42 -20 6];
cfg.atlas            = atlas; % works with atlas!!!
%cfg.funcolorlim     = [-30 30];
raweffect_int.parcel = dummy1.tissue;
ft_sourceplot(cfg,raweffect_int);
raweffect_int.parcel = dummy2.tissue;
ft_sourceplot(cfg,raweffect_int);

% plot averages over Heschl and Superior T. Gyrus
%------------------------------------------------
figure('Position', get(0, 'Screensize'))
idxp(1) = find(contains(parcel.label,{'Heschl_L'}));
idxp(2) = find(contains(parcel.label,{'Heschl_R'}));
idxp(3) = find(contains(parcel.label,{'Temporal_Sup_L'}));
idxp(4) = find(contains(parcel.label,{'Temporal_Sup_R'}));
subplot(2,2,1)
plot(parcel.time*1000,parcel.mom(idxp(1),:),'LineWidth',2,'Color','r');
hold on
plot(parcel.time*1000,parcel.mom(idxp(3),:),'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]);
xlabel('t / ms')
ylabel('mom / a.u.')
grid on
legend('Heschl_L','TemporalSup_L')
subplot(2,2,2)
plot(parcel.time*1000,parcel.mom(idxp(2),:),'LineWidth',2,'Color','b');
hold on
plot(parcel.time*1000,parcel.mom(idxp(4),:),'LineWidth',2,'Color',[0.3010, 0.7450, 0.9330]);
xlabel('t / ms')
ylabel('mom / a.u.')
grid on
legend('Heschl_R','TemporalSup_R')
subplot(2,2,[3,4])
plot(parcel.time*1000,parcel.mom(idxp(1),:),'LineWidth',2,'Color','r');
hold on
plot(parcel.time*1000,parcel.mom(idxp(2),:),'LineWidth',2,'Color','b');
plot(parcel.time*1000,parcel.mom(idxp(3),:),'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]);
plot(parcel.time*1000,parcel.mom(idxp(4),:),'LineWidth',2,'Color',[0.3010, 0.7450, 0.9330]);

xlabel('t / ms')
ylabel('mom / a.u.')
grid on
legend('Heschl_L','Heschl_R','TemporalSup_L','TemporalSup_R')

sgtitle('beamformer gridmodel averaged parcelled data')
set(findall(gcf,'-property','FontSize'),'FontSize',25)
saveas(gcf,fullfile(new_dir,['grandavg_raweffect_heschl',add,add2,'.png']))

%-------------- interpolate whole timeseries on surface-------------------- 
cfg                        = [];
cfg.voxelcoord             = 'no';
cfg.parameter              = 'mom';
cfg.interpmethod           = 'nearest';
raweffect_surfint          = ft_sourceinterpolate(cfg,raweffect,template_surf);
raweffect_surfint.coordsys = 'mni';
save(fullfile(new_dir,['grandavg_raweffect_surfint',add,add2,'.mat']),'raweffect_surfint');

%raweffect_surfint = importdata(fullfile(new_dir,['grandavg_raweffect_surfint',add,add2,'.mat']));
cfg              = [];
% cfg.funparameter = 'mom';
% ft_sourcemovie(cfg,raweffect_surfint);
cfg.parameter = 'mom';
ft_sourceplot_interactive(cfg,raweffect_surfint);
%savefig(gcf,fullfile(new_dir,['grandavg_raweffect_surface',add,add2,'.fig']))

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.latency      = timewin(1,:);
cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,raweffect_surfint);
savefig(gcf,fullfile(new_dir,['grandavg_raweffect_surface_timewin1',add,add2,'.fig']))

%% statistic
cfg                  = [];
cfg.dim              = avg{1}.dim;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.parameter        = 'mom';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.numrandomization = 500;
cfg.alpha            = 0.025; % note that this only implies single-sided testing
cfg.tail             = 0;
cfg.clustertail      = 0;

cfg.design(1,:) = [1:N_subj 1:N_subj];
cfg.design(2,:) = [ones(1,N_subj)*1 ones(1,N_subj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat = ft_sourcestatistics(cfg,avg{:},avg_shuffle{:});

%%
% plot the t values on the MRI 
%-----------------------------
cfg              = [];
cfg.parameter    = 'stat';
statint          = ft_sourceinterpolate(cfg,stat,template_mri);
cfg.parameter    = 'mask';
maskint          = ft_sourceinterpolate(cfg,stat,template_mri);
statint.mask     = maskint.mask;
statint.coordsys = 'mni';

cfg               = [];
% cfg.method      = 'slice';
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.maskparameter = 'mask';
% cfg.funcolorlim = [-5 5];
cfg.funcolormap   = 'jet';
cfg.location      = [64 -32 8];
cfg.atlas         = atlas;
% cfg.latency       = timewin(1,:);
% cfg.avgovertime   = 'yes';
ft_sourceplot(cfg,statint);

%-------------- interpolate whole timeseries on surface-------------------- 
cfg                   = [];
cfg.voxelcoord        = 'no';
cfg.parameter         = {'mask','stat'};
cfg.interpmethod      = 'nearest';
stat_surfint          = ft_sourceinterpolate(cfg,stat,template_surf);
stat_surfint.coordsys = 'mni';


cfg              = [];
% cfg.funparameter = 'stat';
% ft_sourcemovie(cfg,stat_surfint);
cfg.parameter = 'mask';
ft_sourceplot_interactive(cfg,stat_surfint);
cfg.parameter = 'stat';
ft_sourceplot_interactive(cfg,stat_surfint);
% savefig(gcf,fullfile(new_dir,['stat_surface',add,add2,'.fig']))


% my own statistic: signifant clusters and t-values combined 
%-----------------------------------------------------------
my_stat_surfint                    = stat_surfint;
my_stat_surfint.stat               = abs(my_stat_surfint.stat);
my_stat_surfint.significant_points = my_stat_surfint.mask.*my_stat_surfint.stat;

save(fullfile(new_dir,['stat_surfint',add,add2,'.mat']),'my_stat_surfint');
% my_stat_surfint  = importdata(fullfile(new_dir,['stat_surfint',add,add2,'.mat']));

cfg                                = [];
cfg.parameter                      = 'significant_points';
cfg.atlas                          = atlas;
ft_sourceplot_interactive(cfg,my_stat_surfint);
% cfg              = [];
% cfg.funparameter = 'significant_points';
% ft_sourcemovie(cfg,my_stat_surfint);


cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'stat';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.latency      = timewin(1,:);
cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,stat_surfint);
savefig(gcf,fullfile(new_dir,['stat_surface_timewin1',add,add2,'.fig']))

%% Clean up
rmpath(fullfile('Z:','analysis','subject_files'))



% %--------------------------- with parcellation-----------------------------
% % plot the t values on the MRI 
% %-----------------------------
% cfg              = [];
% cfg.parameter    = 'stat';
% statint          = ft_sourceinterpolate(cfg,stat,template_mri_2);
% cfg.parameter    = 'mask';
% maskint          = ft_sourceinterpolate(cfg,stat,template_mri_2);
% statint.mask     = maskint.mask;
% statint.coordsys = 'mni';
% 
% cfg        = [];
% parcel2    = ft_sourceparcellate(cfg,statint, atlas);
% parcelmask = ft_sourceparcellate(cfg,maskint, atlas);
% 
% % create dummy struct
% dummy3    = atlas;
% dummymask = atlas;
% idx1   = dsearchn(raweffect_int_2.time',timewin(1,:)');
% for i=1:size(parcel2.stat,1)
%       dummy3.tissue(find(dummy3.tissue==i))       = mean(parcel.mom(i,idx1(1):idx1(2))); 
%       dummymask.tissue(find(dummymask.tissue==i)) = round(mean(parcelmask.mask(i,idx1(1):idx1(2)))); 
% end
% 
% statint.parcel    = dummy3.tissue;
% statint.coordsys  = 'mni';
% statint.mask      = dummymask.tissue;
% cfg               = [];
% cfg.method        = 'slice';
% cfg.funparameter  = 'parcel';
% cfg.funcolormap   = 'jet';
% cfg.maskparameter = 'mask';
% cfg.renderer      = 'zbuffer';
% % cfg.funcolorlim = [-5 5];
% cfg.atlas         = atlas;
% ft_sourceplot(cfg,statint);