close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots a lot of results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
for i = 3:24
    if i<10; subject='subject0'; else subject='subject'; end
    subjects{i-2} = [subject,num2str(i)]; 
end

% apply ica
ica_status = 1;

% time window for average
timewin = [0.05,0.15];

check = 0;
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'));
% addpath for preprocessing function
addpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));

if ica_status
    add = '_ica';
else 
    add = '';
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
template_surf_infl = importdata(fullfile('Z:','analysis','generated_data','single_subj_T1_1mm','sourcemodel',...
                                         'single_subj_T1_1mm_inflated_4k.mat'));

new_dir = fullfile('Z:','analysis','generated_data','grouplevel','chirps','sourcelevel');
if ~exist(new_dir, 'dir')
   mkdir(new_dir)
end

% load artifical template information for mne
%--------------------------------------------
% eval(subjects{1})
% inflated = importdata(fullfile(subjectdata.sourcemodel,[subjectdata.subjectname,'_inflated_4k.mat']));

%% load data
%--------------------------------------------------------------------------
N_subj                 = length(subjects);
sources_beamf_grid     = cell(1,N_subj);
sources_beamf_grid_raw = cell(1,N_subj);
sources_beamf_surf     = cell(1,N_subj);
sources_mne            = cell(1,N_subj);
dipfit_sources         = cell(1,N_subj);
dipfit_timeseries      = cell(1,N_subj);
dipfit_pos             = zeros(2,3,N_subj);

% load beamformer
%----------------
for s = 1:N_subj
    eval(subjects{s})
    % grid
    %-----
    source_beamf_grid         = importdata(fullfile(subjectdata.chirps_beamformer,...
                                [subjectdata.subjectname,'_beamformer_gridmodel',add,'.mat']));
    sources_beamf_grid{s}     = source_beamf_grid.source;
    sources_beamf_grid{s}.pos = template_grid.pos;
    % generate matrix for moments!!!
    I       = length(sources_beamf_grid{s}.inside);
    moments = NaN(I,length(sources_beamf_grid{s}.time));
    for i = find(sources_beamf_grid{s}.inside)'
            moments(i,:) = sources_beamf_grid{s}.avg.mom{i};   
    end
    sources_beamf_grid_raw{s}         = sources_beamf_grid{s};
    sources_beamf_grid_raw{s}.avg.mom = moments;
    sources_beamf_grid{s}.avg.mom     = abs(moments);
    % normalize
    sources_beamf_grid_raw{s}.avg.mom = sources_beamf_grid_raw{s}.avg.mom./max(abs(sources_beamf_grid_raw{s}.avg.mom),[],'all');
    sources_beamf_grid{s}.avg.mom     = sources_beamf_grid{s}.avg.mom./max(abs(sources_beamf_grid{s}.avg.mom),[],'all');
    
    % surface
    %--------
    source_beamf_surf             = importdata(fullfile(subjectdata.chirps_beamformer,...
                                    [subjectdata.subjectname,'_beamformer_surfacemodel',add,'.mat']));
    sources_beamf_surf{s}         = source_beamf_surf.source;
    sources_beamf_surf{s}.pos     = template_surf.pos; 
    % generate matrix for moments!!!
    I       = length(sources_beamf_surf{s}.inside);
    moments = NaN(I,length(sources_beamf_surf{s}.time));
    for i = find(sources_beamf_surf{s}.inside)'
            moments(i,:) = sources_beamf_surf{s}.avg.mom{i};   
    end
    sources_beamf_surf{s}.avg.mom = abs(moments);
    % normalize
    sources_beamf_surf{s}.avg.mom = sources_beamf_surf{s}.avg.mom./max(abs(sources_beamf_surf{s}.avg.mom),[],'all');   
end

% mne
%----
% projected moments seems to be less focal -> use norm of moments
for s = 1:N_subj
    eval(subjects{s})
    source_mne             = importdata(fullfile(subjectdata.chirps_mne,[subjectdata.subjectname,'_mne_sources',add,'.mat']));
    % project the source to its strongest orientation, i.e. the direction that explains most of the source variance. 
    % That is equivalent to taking the largest eigenvector of the source timeseries.
%     cfg                    = [];
%     cfg.projectmom         = 'yes';
%     sources_mne{s}         = ft_sourcedescriptives(cfg,source_mne.source);  
    sources_mne{s}         = source_mne.source;
    sources_mne{s}.pos     = template_surf.pos;
    % generate matrix for moments!!!
    I       = length(sources_mne{s}.inside);
    moments = NaN(I,length(sources_mne{s}.time));
    for i = find(sources_mne{s}.inside)'
            %moments(i,:) = sources_mne{s}.avg.mom{i};  
            moments(i,:) = sqrt(sum(sources_mne{s}.avg.mom{i}.^2,1)); 
    end
    sources_mne{s}.avg.mom = abs(moments);
    % normalize
    sources_mne{s}.avg.mom = sources_mne{s}.avg.mom./max(abs(sources_mne{s}.avg.mom),[],'all');
end

% dipolefitting
%--------------
for s = 1:N_subj
    eval(subjects{s})
    source_dipfit = importdata(fullfile(subjectdata.chirps_dipolefitting,[subjectdata.subjectname,...
                    '_results_dipolefitting',add,'.mat']));
    dipfit_pos(:,:,s)         = source_dipfit.pos_nosym_mni;
    dipfit_sources{s}         = source_dipfit.source_nosym;
    dipfit_sources{s}.dip.pos = source_dipfit.pos_nosym_mni;
    dipfit_timeseries{s}      = source_dipfit.source_timeseries;
end

%% calculate grand average  
%--------------------------------------------------------------------------

% beamformer
%-----------
cfg                     = [];
cfg.latency             = 'all';
cfg.parameter           = 'mom';
grandavg_beamf_grid     = ft_sourcegrandaverage(cfg,sources_beamf_grid{:});
grandavg_beamf_grid_raw = ft_sourcegrandaverage(cfg,sources_beamf_grid_raw{:});

cfg                     = [];
cfg.latency             = 'all';
cfg.parameter           = 'mom';
grandavg_beamf_surf     = ft_sourcegrandaverage(cfg,sources_beamf_surf{:});
grandavg_beamf_surf.tri = template_surf_infl.tri;
grandavg_beamf_surf.pos = template_surf_infl.pos;

% mne
%----
cfg              = [];
cfg.latency      = 'all';
cfg.parameter    = 'mom';
grandavg_mne     = ft_sourcegrandaverage(cfg, sources_mne{:});
grandavg_mne.tri = template_surf_infl.tri;
grandavg_mne.pos = template_surf_infl.pos;

%% Plot data
%--------------------------------------------------------------------------

%%% beamformer - grid - raw moments %%%
%--------------------------------------
% interpolate onto mri
cfg                                  = [];
cfg.voxelcoord                       = 'no';
cfg.parameter                        = 'mom';
cfg.interpmethod                     = 'nearest';
grandavg_beamf_grid_raw_int          = ft_sourceinterpolate(cfg,grandavg_beamf_grid_raw,template_mri);
grandavg_beamf_grid_raw_int.coordsys = 'mni';

% Plot results
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.atlas        = atlas;
cfg.latency      = timewin;
cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,grandavg_beamf_grid_raw_int);

%-----------------------grid with parcellation-----------------------------
cfg                           = [];
cfg.voxelcoord                = 'no';
cfg.parameter                 = 'mom';
cfg.interpmethod              = 'nearest';
grandavg_beamf_grid_raw_int_2 = ft_sourceinterpolate(cfg, grandavg_beamf_grid_raw, template_mri_2);
grandavg_beamf_grid_raw_int_2.coordsys = 'mni';

cfg     = [];
parcel2 = ft_sourceparcellate(cfg,grandavg_beamf_grid_raw_int_2,atlas);

% create a dummy structure where we identify the moment values per voxel
% but remember - moments is a timeseries
dummy = atlas;
idx   = dsearchn(grandavg_beamf_grid_raw_int_2.time',timewin');
for i=1:size(parcel2.mom,1)
      dummy.tissue(find(dummy.tissue==i)) = mean(parcel2.mom(i,idx(1):idx(2))); 
end

grandavg_beamf_grid_raw_int_2.parcel   = dummy.tissue;
grandavg_beamf_grid_raw_int_2.coordsys = 'mni';
cfg                                    = [];
cfg.method                             = 'ortho';
cfg.funparameter                       = 'parcel';
cfg.funcolormap                        = 'jet';
cfg.renderer                           = 'zbuffer';
cfg.location                           = [-42 -20 6];
cfg.atlas                              = atlas; % works with atlas!!!
%cfg.funcolorlim                       = [-30 30];
ft_sourceplot(cfg,grandavg_beamf_grid_raw_int_2);

% plot averages over Heschl
figure
idxp(1) = find(contains(parcel2.label,{'Heschl_L'}));
idxp(2) = find(contains(parcel2.label,{'Heschl_R'}));
subplot(2,2,1)
plot(parcel2.time*1000,parcel2.mom(idxp(1),:),'LineWidth',2,'Color','r');
xlabel('t / ms')
ylabel('mom / a.u.')
grid on
legend('Heschl_L')
subplot(2,2,2)
plot(parcel2.time*1000,parcel2.mom(idxp(2),:),'LineWidth',2,'Color','b');
xlabel('t / ms')
ylabel('mom / a.u.')
grid on
legend('Heschl_R')
subplot(2,2,[3,4])
plot(parcel2.time*1000,parcel2.mom(idxp(1),:),'LineWidth',2,'Color','r');
hold on
plot(parcel2.time*1000,parcel2.mom(idxp(2),:),'LineWidth',2,'Color','b');
xlabel('t / ms')
ylabel('mom / a.u.')
grid on
legend('Heschl_L','Heschl_R')
sgtitle('beamformer gridmodel averaged parcelled data')
set(findall(gcf,'-property','FontSize'),'FontSize',25)
%--------------------------------------------------------------------------

%%% beamformer - grid %%%
%------------------------
% interpolate onto mri
cfg                              = [];
cfg.voxelcoord                   = 'no';
cfg.parameter                    = 'mom';
cfg.interpmethod                 = 'nearest';
grandavg_beamf_grid_int          = ft_sourceinterpolate(cfg,grandavg_beamf_grid,template_mri);
grandavg_beamf_grid_int.coordsys = 'mni';

% Plot results
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.atlas        = atlas;
cfg.latency      = timewin;
cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,grandavg_beamf_grid_int);

%-----------------------grid with parcellation-----------------------------
cfg                       = [];
cfg.voxelcoord            = 'no';
cfg.parameter             = 'mom';
cfg.interpmethod          = 'nearest';
grandavg_beamf_grid_int_2 = ft_sourceinterpolate(cfg, grandavg_beamf_grid, template_mri_2);
grandavg_beamf_grid_int_2.coordsys = 'mni';

cfg    = [];
parcel = ft_sourceparcellate(cfg,grandavg_beamf_grid_int_2,atlas);

% create a dummy structure where we identify the moment values per voxel
% but remember - moments is a timeseries
dummy = atlas;
idx   = dsearchn(grandavg_beamf_grid_int_2.time',timewin');
for i=1:size(parcel.mom,1)
      dummy.tissue(find(dummy.tissue==i)) = mean(parcel.mom(i,idx(1):idx(2))); 
end

grandavg_beamf_grid_int_2.parcel   = dummy.tissue;
grandavg_beamf_grid_int_2.coordsys = 'mni';
cfg                                = [];
cfg.method                         = 'ortho';
cfg.funparameter                   = 'parcel';
cfg.funcolormap                    = 'jet';
cfg.renderer                       = 'zbuffer';
cfg.location                       = [-42 -20 6];
cfg.atlas                          = atlas; % works with atlas!!!
%cfg.funcolorlim                   = [-30 30];
ft_sourceplot(cfg,grandavg_beamf_grid_int_2);
savefig(gcf,fullfile(new_dir,['grandavg_beamf_grid_parcelled',add,'.fig']))

% plot averages over Heschl
figure
idxp(1) = find(contains(parcel.label,{'Heschl_L'}));
idxp(2) = find(contains(parcel.label,{'Heschl_R'}));
subplot(2,2,1)
plot(parcel.time*1000,parcel.mom(idxp(1),:),'LineWidth',2,'Color','r');
xlabel('t / ms')
ylabel('|mom| / a.u.')
grid on
legend('Heschl_L')
subplot(2,2,2)
plot(parcel.time*1000,parcel.mom(idxp(2),:),'LineWidth',2,'Color','b');
xlabel('t / ms')
ylabel('|mom| / a.u.')
grid on
legend('Heschl_R')
subplot(2,2,[3,4])
plot(parcel.time*1000,parcel.mom(idxp(1),:),'LineWidth',2,'Color','r');
hold on
plot(parcel.time*1000,parcel.mom(idxp(2),:),'LineWidth',2,'Color','b');
xlabel('t / ms')
ylabel('|mom| / a.u.')
grid on
legend('Heschl_L','Heschl_R')
sgtitle('beamformer gridmodel averaged parcelled data')
set(findall(gcf,'-property','FontSize'),'FontSize',25)
%--------------------------------------------------------------------------

% template_surface = ft_read_headshape(['Z:\Software\MEG_EEG_Toolboxen\fieldtrip-20191127\' ... 
%                             'template\sourcemodel\cortex_5124.surf.gii']);
%-------------- interpolate whole timeseries on surface-------------------- 
cfg                                  = [];
cfg.voxelcoord                       = 'no';
cfg.parameter                        = 'mom';
cfg.interpmethod                     = 'nearest';
grandavg_beamf_grid_surfint          = ft_sourceinterpolate(cfg,grandavg_beamf_grid,template_surf);
grandavg_beamf_grid_surfint.coordsys = 'mni';

cfg              = [];
cfg.funparameter = 'mom';
ft_sourcemovie(cfg,grandavg_beamf_grid_surfint);
% cfg.parameter = 'mom';
% ft_sourceplot_interactive(cfg,grandavg_beamf_grid_surfint);

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.latency      = timewin;
cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,grandavg_beamf_grid_surfint);
savefig(gcf,fullfile(new_dir,['grandavg_beamf_grid_intonsurface',add,'.fig']))

%--------------- project averaged data on brain surface--------------------
grandavg_beamf_grid_int_surf     = grandavg_beamf_grid_int;
grandavg_beamf_grid_int_surf.mom = mean(grandavg_beamf_grid_int_surf.mom(:,idx(1):idx(2)),2);

cfg              = [];
cfg.method       = 'surface';
cfg.atlas        = atlas;
cfg.funparameter = 'mom';
%cfg.funcolorlim = [-30 30];
cfg.funcolormap  = 'jet';
cfg.projmethod   = 'nearest';
cfg.surfinflated = 'surface_inflated_both_caret.mat';
cfg.projthresh   = 0.7; % in percent - adapt threshold!!!
cfg.camlight     = 'no';
ft_sourceplot(cfg, grandavg_beamf_grid_int_surf);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull

%%% beamformer - surface %%%
%---------------------------
cfg              = [];
cfg.funparameter = 'mom';
ft_sourcemovie(cfg,grandavg_beamf_surf);

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.latency      = timewin;
cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,grandavg_beamf_surf);
savefig(gcf,fullfile(new_dir,['grandavg_beamf_surf',add,'.fig']))

%%% mne %%%
%----------
cfg              = [];
cfg.atlas        = atlas;
% cfg.funparameter = 'mom';
% ft_sourcemovie(cfg,grandavg_mne);
cfg.parameter = 'mom';
ft_sourceplot_interactive(cfg,grandavg_mne); % works with atlas!
savefig(gcf,fullfile(new_dir,['grandavg_mne_movie',add,'.fig']))

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'mom';
cfg.location     = [64 -32 8];
cfg.funcolormap  = 'jet';
cfg.latency      = timewin;
cfg.avgovertime  = 'yes';
ft_sourceplot(cfg,grandavg_mne);
savefig(gcf,fullfile(new_dir,['grandavg_mne',add,'.fig']))

%%
%%% dipolefitting %%%
%--------------------
[bad,label_left,label_right] = check_dipolefitting(dipfit_pos,subjects);
[meani,meantimecourse]       = dipolefitting_mean_values(dipfit_sources,dipfit_timeseries,bad,label_left);

% plot 3D dipoles
color = parula(N_subj);
figure
hold on
for s = setdiff([1:N_subj],bad)
ft_plot_dipole(dipfit_sources{s}.dip.pos(1,:), mean(dipfit_sources{s}.dip.mom(1:3,:),2), 'color', color(s,:),  'unit', 'mm','alpha',0.3)
ft_plot_dipole(dipfit_sources{s}.dip.pos(2,:), mean(dipfit_sources{s}.dip.mom(4:6,:),2), 'color', color(s,:),  'unit', 'mm','alpha',0.3)
end
ft_plot_dipole(meani.dip.pos(1,:), mean(dipfit_sources{s}.dip.mom(1:3,:),2), 'color',[0.6350, 0.0780, 0.1840], 'unit', 'mm','thickness',5)
ft_plot_dipole(meani.dip.pos(2,:), mean(dipfit_sources{s}.dip.mom(4:6,:),2), 'color', [0.6350, 0.0780, 0.1840], 'unit', 'mm','thickness',5)

pos = mean(meani.dip.pos,1);
ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)
savefig(gcf,fullfile(new_dir,['grandavg_dipolefitting',add,'.fig']))

pos_dip_l = my_atlas_lookup(atlas,meani.dip.pos(1,:),'coordsys','mni');
pos_dip_r = my_atlas_lookup(atlas,meani.dip.pos(2,:),'coordsys','mni');
disp('Dipolefitting')
disp('-------------')
disp(['mean dipole position left: ',pos_dip_l{1}])
disp(['mean dipole position right: ',pos_dip_r{1}])

% plot averaged timecourses
%--------------------------
figure
c = {'b','g','r'};
subplot(2,1,1); 
title('left hemisphere')
hold on
arrayfun(@(i) plot(dipfit_timeseries{1}.time*1000,meantimecourse{1}(i,:),'-','color',c{i},'LineWidth',2),1:3)
p1          = plot(dipfit_timeseries{1}.time*1000,sqrt(sum(meantimecourse{1}.^2,1)),'k-','LineWidth',3);
p1.Color(4) = 0.25;
xlabel('t / ms')
ylabel('mom / a.u.')
legend({'x', 'y', 'z','|.|'});
grid on
subplot(2,1,2); 
title('right hemisphere')
hold on
arrayfun(@(i) plot(dipfit_timeseries{2}.time,meantimecourse{2}(i,:),'-','color',c{i},'LineWidth',2),1:3)
p1          = plot(dipfit_timeseries{1}.time,sqrt(sum(meantimecourse{2}.^2,1)),'k-','LineWidth',3);
p1.Color(4) = 0.25;
xlabel('t / ms')
ylabel('mom / a.u.')
legend({'x', 'y', 'z','|.|'});
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',25)

%% Plot dipole timecourse together with beamforming gridmodel timecourse
figure('Position', get(0, 'Screensize'))
%figure
c = {'b','g','r'};
subplot(2,2,1); 
title('left hemisphere')
hold on
arrayfun(@(i) plot(dipfit_timeseries{1}.time*1000,meantimecourse{1}(i,:),'-','color',c{i},'LineWidth',2),1:3)
p1          = plot(dipfit_timeseries{1}.time*1000,sqrt(sum(meantimecourse{1}.^2,1)),'k-','LineWidth',3);
p1.Color(4) = 0.25;
xlabel('t / ms')
ylabel('mom / a.u.')
xlim([-250,500])
legend({'x', 'y', 'z','|.|'});
grid on
h = text(-400,-0.5,'dipole fitting','FontSize',25,'FontWeight','bold');
set(h,'Rotation',90);

subplot(2,2,2); 
title('right hemisphere')
hold on
arrayfun(@(i) plot(dipfit_timeseries{2}.time*1000,meantimecourse{2}(i,:),'-','color',c{i},'LineWidth',2),1:3)
p1          = plot(dipfit_timeseries{1}.time*1000,sqrt(sum(meantimecourse{2}.^2,1)),'k-','LineWidth',3);
p1.Color(4) = 0.25;
xlabel('t / ms')
ylabel('mom / a.u.')
xlim([-250,500])
legend({'x', 'y', 'z','|.|'});
grid on
% plot averages over Heschl
idxp(1) = find(contains(parcel2.label,{'Heschl_L'}));
idxp(2) = find(contains(parcel2.label,{'Heschl_R'}));
subplot(2,2,3)
plot(parcel2.time*1000,parcel2.mom(idxp(1),:),'LineWidth',2,'Color','r');
xlabel('t / ms')
ylabel('mom / a.u.')
xlim([-250,500])
grid on
legend('Heschl_L')
h = text(-400,-0.03,'beamforming','FontSize',25,'FontWeight','bold');
set(h,'Rotation',90);

subplot(2,2,4)
plot(parcel2.time*1000,parcel2.mom(idxp(2),:),'LineWidth',2,'Color','b');
xlabel('t / ms')
ylabel('mom / a.u.')
xlim([-250,500])
grid on
legend('Heschl_R')
sgtitle('source waveforms')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
saveas(gcf,fullfile(new_dir,['grandavg_dipolefitting_beamforming_timecourse',add,'.png']))

%% old
% for s = 1:N_subj
%     eval(subjects{s})
%     % beamformer
%     %-----------
%     source_beamf              = importdata(fullfile(subjectdata.chirps_beamformer,...
%                                 [subjectdata.subjectname,'_beamformer_gridmodel.mat']));
%   % sources_beamf{s}          = source_beamf.source;
%     cfg                       = [];
%     cfg.operation             = 'x1.^2';
%     cfg.parameter             = 'mom';
%     sourcepre_beamf           = ft_math(cfg,source_beamf.sourcepre);
%     sourcepst_beamf           = ft_math(cfg,source_beamf.sourcepst);
%     
%     cfg                       = [];
%     cfg.parameter             = 'mom';
%     cfg.operation             = 'subtract';
%     sources_diff_beamf{s}     = ft_math(cfg,sourcepst_beamf,sourcepre_beamf);
%     
%     % generate matrix for moments!!!
%     I       = length(sources_diff_beamf{s}.inside);
%     moments = NaN(I,length(sources_diff_beamf{s}.time));
%     for i = find(sources_diff_beamf{s}.inside)'
%             moments(i,:) = sources_diff_beamf{s}.mom{i};   
%     end
%     sources_diff_beamf{s}.mom = moments;
%     sources_diff_beamf{s}.mom = sources_diff_beamf{s}.mom./max(abs(sources_diff_beamf{s}.mom),[],'all');
%     
%     sources_diff_beamf{s}.pos = template_grid.pos;
%     info_beamf{s}             = source_beamf.info;  
%     kappas_beamf(s,1:2)       = info_beamf{s}.kappa_noise;
%     kappas_beamf(s,3:4)       = info_beamf{s}.kappa_data;
%     
%     % mne
%     %----
%     source_mne             = importdata(fullfile(subjectdata.chirps_mne,[subjectdata.subjectname,'_mne_sources.mat']));
%     % project the source to its strongest orientation, i.e. the direction that explains most of the source variance. 
%     % That is equivalent to taking the largest eigenvector of the source timeseries.
%     cfg                    = [];
%     cfg.projectmom         = 'yes';
%     sources_mne{s}         = ft_sourcedescriptives(cfg,source_mne.source);  
%     sources_mne{s}.avg.pow = sources_mne{s}.avg.pow./max(abs(sources_mne{s}.avg.pow),[],'all');
%     sources_mne{s}.pos     = sources_mne{1}.pos;
%     info_mne{s}            = source_mne.info;  
%     kappas_mne(s,1:2)      = info_mne{s}.kappa_noise;
%     kappas_mne(s,3:4)      = info_mne{s}.kappa_data;
%  
% end
% kappas_beamf = array2table(kappas_beamf,'VariableNames',{'noise-mag','noise-grad','data-mag','data-grad'});
% kappas_mne   = array2table(kappas_mne,'VariableNames',{'noise-mag','noise-grad','data-mag','data-grad'});




%% Clean up
rmpath(fullfile('Z:','analysis','subject_files'));
rmpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));

%% functions
function [bad_subjects,label_left,label_right] = check_dipolefitting(dipfit_pos,subjects)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% searches for bad dipole fits and returns bad subject numbers and
% labelling of left and right data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dipfit_label = subjects;
N_subj       = length(subjects);
dipfit_bad   = {};
label_left   = zeros(1,N_subj);
label_right  = zeros(1,N_subj);
% 1.) remove subjects with dipoles in same hemissphere 
% -> look if x-coordinates has different signs
xcoord               = squeeze(dipfit_pos(:,1,:));
idx1                 = find(prod(xcoord,1) > 0); 
dipfit_pos(:,:,idx1) = [];

if ~isempty(idx1)
    dipfit_bad         = [dipfit_bad,dipfit_label(idx1)];
    dipfit_label(idx1) = [];
    label_left(idx1)   = [];
    label_right(idx1)  = [];
end
% 2.) sort positions in left and right
idx2         = squeeze(dipfit_pos(:,1,:) < 0);
I            = length(dipfit_label);
dipfit_pos_l = zeros(I,3);
dipfit_pos_r = zeros(I,3);

for i = 1:I
    idx_l = find(idx2(:,i));  % <0
    idx_r = find(~idx2(:,i)); % >0
    
    label_left(i)   = idx_l;
    label_right(i)  = idx_r;
    
    dipfit_pos_l(i,:) = dipfit_pos(idx_l,:,i);
    dipfit_pos_r(i,:) = dipfit_pos(idx_r,:,i);
end
% 3.) remove outliers with distanc metric and MAD
distance_l = zeros(I);
distance_r = zeros(I);
for i = 1:I
    for j = 1:I
        distance_l(i,j) = norm(dipfit_pos_l(i,:)-dipfit_pos_l(j,:));
        distance_r(i,j) = norm(dipfit_pos_r(i,:)-dipfit_pos_r(j,:));
    end
end
% new scalar distance metric
distance_sum_l = sum(distance_l,1);
distance_sum_r = sum(distance_r,1);

% by default, an outlier is a value that is more than three scaled 
% median absolute deviations (MAD) away from the median
idx_l    = find(isoutlier(distance_sum_l)); 
idx_r    = find(isoutlier(distance_sum_r)); 
combined = union(idx_l,idx_r);
if ~isempty(combined)
        dipfit_bad             = [dipfit_bad,dipfit_label(combined)];
        dipfit_label(combined) = [];
        label_left(combined)   = [];
        label_right(combined)  = [];
end

% 4.) get numbers of bad subjects
I = length(dipfit_bad);
for i = 1:I
    bad_subjects = find(contains(subjects,dipfit_bad));
end

end

%--------------------------------------------------------------------------

function [meani,meantimecourse] = dipolefitting_mean_values(dipfit_sources,dipfit_timeseries,bad,label_left)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generates mean values from good subjects for dipole fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanpos        = cell(1,2);
meanmom        = cell(1,2);
meantimecourse = cell(1,2);
N_subj         = length(dipfit_sources);
pos_l          = zeros(N_subj,3);
pos_r          = zeros(N_subj,3);
n              = 1;

% only loop over good subjects
for s = setdiff([1:N_subj],bad) 
    if label_left(n)==1 
        % mean moments and positions for each hemisphere
        pos_l(n,:) = dipfit_sources{s}.dip.pos(1,:);
        pos_r(n,:) = dipfit_sources{s}.dip.pos(2,:);
        
        dipfit_mom_l(:,:,n) = dipfit_sources{s}.dip.mom(1:3,:);
        dipfit_mom_r(:,:,n) = dipfit_sources{s}.dip.mom(4:6,:);
        
        % mean timecourses
        dipfit_timecourse_l(:,:,n) = dipfit_timeseries{s}.dip.mom(1:3,:);
        dipfit_timecourse_r(:,:,n) = dipfit_timeseries{s}.dip.mom(4:6,:);      
    else
        % mean moments and positions for each hemisphere
        pos_l(n,:) = dipfit_sources{s}.dip.pos(2,:);
        pos_r(n,:)  = dipfit_sources{s}.dip.pos(1,:);
        
        dipfit_mom_l(:,:,n) = dipfit_sources{s}.dip.mom(1:3,:);
        dipfit_mom_r(:,:,n) = dipfit_sources{s}.dip.mom(4:6,:);
        
        % mean timecourses
        dipfit_timecourse_l(:,:,n) = dipfit_timeseries{s}.dip.mom(1:3,:);
        dipfit_timecourse_r(:,:,n) = dipfit_timeseries{s}.dip.mom(4:6,:);
    end  
    % optional: normalize
    dipfit_timecourse_l(:,:,n) = dipfit_timecourse_l(:,:,n)./max(abs(dipfit_timecourse_l(:,:,n)),[],'all');
    dipfit_timecourse_r(:,:,n) = dipfit_timecourse_r(:,:,n)./max(abs(dipfit_timecourse_r(:,:,n)),[],'all');
        
    n = n+1;
end

% 1. left, 2. right
meanpos{1} = mean(pos_l,1);
meanpos{2} = mean(pos_r,1);

meanmom{1} = mean(dipfit_mom_l,3);
meanmom{2} = mean(dipfit_mom_r,3);

meantimecourse{1} = mean(dipfit_timecourse_l,3);
meantimecourse{2} = mean(dipfit_timecourse_r,3);

meani         = dipfit_sources{1};
meani.dip.pos = [meanpos{1};meanpos{2}];
meani.dip.mom = [meanmom{1};meanmom{2}];

end