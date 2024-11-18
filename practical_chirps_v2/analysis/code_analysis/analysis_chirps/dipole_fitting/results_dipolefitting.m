close all; clear all; clc;

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject04'};

% choose files
files2preproc = 'stories_maxfilter';
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(['Z:' filesep 'analysis' filesep 'subject_files']);

%% initialize

% loop over subjects 
%-------------------
N_subj = length(subjects);

for s = 1:N_subj
    
% subject selected
subject = subjects{s};
eval(subject)

headmodel_meg = importdata([subjectdata.headmodel filesep 'headmodel_meg.mat']); % mm
mri_segmented = importdata([subjectdata.headmodel filesep 'mri_segmented.mat']); % mm

results = importdata([subjectdata.chirps_dipolefitting filesep files2preproc ...
                     filesep 'results_dipolefitting.mat']);  

%% plot dipole positions
figure
hold on
ft_plot_dipole(results.source_planar.dip.pos(1,:), mean(results.source_planar.dip.mom(1:3,:),2), 'color', 'g', 'unit', 'mm')
ft_plot_dipole(results.source_planar.dip.pos(2,:), mean(results.source_planar.dip.mom(4:6,:),2), 'color', 'g', 'unit', 'mm')

ft_plot_dipole(results.source_planar_nosym.dip.pos(1,:), mean(results.source_planar_nosym.dip.mom(1:3,:),2), 'color', 'm',  'unit', 'mm')
ft_plot_dipole(results.source_planar_nosym.dip.pos(2,:), mean(results.source_planar_nosym.dip.mom(4:6,:),2), 'color', 'm',  'unit', 'mm')

pos = mean(results.source_planar_nosym.dip.pos,1);
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)

%% moving dipole model
figure
for i=1:numel(results.moving_source.dip)
    pos1(i,:) = results.moving_source.dip(i).pos(1,:);
    pos2(i,:) = results.moving_source.dip(i).pos(2,:);
end
hold on
plot3(pos1(:,1), pos1(:,2), pos1(:,3),'r.')
plot3(pos2(:,1), pos2(:,2), pos2(:,3),'g.')
pos = (mean(pos1, 1) + mean(pos2, 1))/2;   

ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off

%% source timecourses

figure
c = {'b','g','r'};
subplot(2,1,1); title('megplanar: probably right')
hold on
arrayfun(@(i) plot(results.source_timeseries.time, results.source_timeseries.dip.mom(i,:),'-','color',c{i}),1:3)
plot(results.source_timeseries.time, sqrt(sum(results.source_timeseries.dip.mom(1:3,:).^2,1)),'k-','LineWidth',2)
legend({'x', 'y', 'z','|.|'});
%ylim([-0.01 0.02])
grid on

subplot(2,1,2); title('megplanar: probably left')
hold on
%plot(source_all.time, source_all.dip.mom(4:6,:), '-')
arrayfun(@(i) plot(results.source_timeseries.time,results.source_timeseries.dip.mom(i,:),'-','color',c{i-3}),4:6)
plot(results.source_timeseries.time, sqrt(sum(results.source_timeseries.dip.mom(4:6,:).^2,1)),'k-','LineWidth',2)
legend({'x', 'y', 'z','|.|'});
%ylim([-0.01 0.02])
grid on


end

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
