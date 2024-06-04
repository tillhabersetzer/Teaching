close all; clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare empty room measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
%--------------------------------------------------------------------------
% choose subject 
subject = 'subject04';
%--------------------------------------------------------------------------

% addpath for subject_files information
addpath(fullfile('Z:','analysis','subject_files'));
% addpath for ica functions
addpath(fullfile('Z:','analysis','preprocessing_batch','helper_functions'));
% addpath for preprocessing function
addpath(fullfile('Z:','analysis','analysis_chirps','helper_functions'));

eval(subject)

%% analysis
filenames  = horzcat(get_filenames(subjectdata,'empty_pre'),get_filenames(subjectdata,'empty_post'));
noise      = give_noise(filenames,subjectdata);

filenames  = horzcat(get_filenames(subjectdata,'empty_pre_maxfilter'),get_filenames(subjectdata,'empty_post_maxfilter'));
noise_tsss = give_noise(filenames,subjectdata);

filenames = strrep(filenames,'tsss','sss');
noise_sss = give_noise(filenames,subjectdata);

cfg                  = [];
cfg.channel          = 'meg';
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg,noise);
avg_noise_tsss       = ft_timelockanalysis(cfg,noise_tsss);
avg_noise_sss        = ft_timelockanalysis(cfg,noise_sss);

% noise covariance matrix
%------------------------
selmag  = ft_chantype(avg_noise.label, 'megmag');
selgrad = ft_chantype(avg_noise.label, 'megplanar');
C = avg_noise.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
figure
subplot(1,3,1)
imagesc(C);
hold on;
plot(102.5.*[1 1],[0 306],'w','linewidth',2);
plot([0 306],102.5.*[1 1],'w','linewidth',2);
title('noise')
C = avg_noise_tsss.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
subplot(1,3,2)
imagesc(C);
hold on;
plot(102.5.*[1 1],[0 306],'w','linewidth',2);
plot([0 306],102.5.*[1 1],'w','linewidth',2);
title('noise tsss')
C = avg_noise_sss.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
subplot(1,3,3)
imagesc(C);
hold on;
plot(102.5.*[1 1],[0 306],'w','linewidth',2);
plot([0 306],102.5.*[1 1],'w','linewidth',2);
title('noise sss')
sgtitle('noise covariance matrix')

%%
kappa      = give_kappa_value(avg_noise.cov,avg_noise.label,{'megmag','megplanar'});
kappa_sss  = give_kappa_value(avg_noise_sss.cov,avg_noise_sss.label,{'megmag','megplanar'});
kappa_tsss = give_kappa_value(avg_noise_tsss.cov,avg_noise_tsss.label,{'megmag','megplanar'});

% singular values
%----------------
[~,s,~] = svd(avg_noise.cov(selgrad,selgrad));
figure
subplot(1,3,1)
plot(log10(diag(s)),'o');
grid on
title('noise')
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');
text(0.9*xlim(2),ylim(2),['mag: ',num2str(kappa(1)),' grad: ',num2str(kappa(2))])
[~,s,~] = svd(avg_noise_tsss.cov);
subplot(1,3,2)
plot(log10(diag(s)),'o');
grid on
title('noise tsss')
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');
text(0.9*xlim(2),ylim(2),['mag: ',num2str(kappa_sss(1)),' grad: ',num2str(kappa_sss(2))])
[~,s,~] = svd(avg_noise_sss.cov);
subplot(1,3,3)
plot(log10(diag(s)),'o');
grid on
title('noise sss')
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');
text(0.9*xlim(2),ylim(2),['mag: ',num2str(kappa_tsss(1)),' grad: ',num2str(kappa_tsss(2))])
sgtitle('singulare values')

%% Clean up
rmpath(['Z:' filesep 'analysis' filesep 'subject_files'])
rmpath(['Z:' filesep 'analysis' filesep 'preprocessing_batch' filesep 'helper_functions'])
rmpath(['Z:' filesep 'analysis' filesep 'analysis_chirps' filesep 'helper_functions']);