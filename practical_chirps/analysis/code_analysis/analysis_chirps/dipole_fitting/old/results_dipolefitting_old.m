close all; clear all; clc;

% Settings
%--------------------------------------------------------------------------
% choose subject 
subjects = {'subject04'};
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

headmodel_meg    = importdata([subjectdata.headmodel filesep 'headmodel_meg.mat']); % mm
mri_segmented    = importdata([subjectdata.headmodel filesep 'mri_segmented.mat']); % mm

% find directories
listing = dir(subjectdata.chirps_dipolefitting);
folders = {listing(3:end).name};
N_files = length(folders);
results = cell(1,N_files);
% colors
c = jet(N_files);

for n = 1:N_files
    results{n} = importdata([subjectdata.chirps_dipolefitting filesep folders{n} ...
                             filesep 'results_dipolefitting.mat']);  
end

%% plot dipole positions
figure
for n = 1:N_files
ft_plot_dipole(results{n}.source_planar_nosym.dip.pos(1,:), mean(results{n}.source_planar_nosym.dip.mom(1:3,:),2), 'color', c(n,:));
ft_plot_dipole(results{n}.source_planar_nosym.dip.pos(2,:), mean(results{n}.source_planar_nosym.dip.mom(4:6,:),2), 'color', c(n,:),'visible','off')
end
posi = mean(results{1}.source_planar_nosym.dip.pos,1);
ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', posi, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', posi, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', posi, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(posi, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)

%% moving dipole model
figure
for n = 1:N_files
    subplot(2,3,n)
    for i=1:numel(results{n}.moving_source.dip)
        pos1(i,:) = results{n}.moving_source.dip(i).pos(1,:);
        pos2(i,:) = results{n}.moving_source.dip(i).pos(2,:);
    end
    hold on
    plot3(pos1(:,1), pos1(:,2), pos1(:,3),'.','color',c(n,:))
    plot3(pos2(:,1), pos2(:,2), pos2(:,3),'*','color',c(n,:))
    %pos = (mean(pos1, 1) + mean(pos2, 1))/2;   

    ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', posi, 'orientation', [1 0 0], 'resolution', 0.1)
    ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', posi, 'orientation', [0 1 0], 'resolution', 0.1)
    ft_plot_slice(mri_segmented_cm.anatomy, 'transform', mri_segmented_cm.transform, 'location', posi, 'orientation', [0 0 1], 'resolution', 0.1)

    ft_plot_crosshair(posi, 'color', [1 1 1]/2);
    axis tight
    axis off
    
    title(folders{n})
end

%% source timecourses

figure
for n = 1:N_files
subplot(N_files,2,2*n-1); title([folders{n} '| left'])
hold on
plot(results{n}.source_timeseries.time, results{n}.source_timeseries.dip.mom(1:3,:), '-')
plot(results{n}.source_timeseries.time, sqrt(sum(results{n}.source_timeseries.dip.mom(1:3,:).^2,1)),'k-','LineWidth',2)
legend({'x', 'y', 'z','|.|'});
%ylim([-0.01 0.02])
grid on

subplot(N_files,2,2*n); title([folders{n} ' | right'])
hold on
plot(results{n}.source_timeseries.time, results{n}.source_timeseries.dip.mom(4:6,:), '-')
plot(results{n}.source_timeseries.time, sqrt(sum(results{n}.source_timeseries.dip.mom(4:6,:).^2,1)),'k-','LineWidth',2)
legend({'x', 'y', 'z','|.|'});
%ylim([-0.01 0.02])
grid on
end

end

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
