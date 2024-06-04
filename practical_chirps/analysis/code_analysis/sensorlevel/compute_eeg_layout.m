%--------------------------------------------------------------------------
% layout
%--------------------------------------------------------------------------

close all
clear 
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','code_analysis','subjectdata'))
eval('main_settings')

addpath(fullfile('..','code_analysis','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
% choose subject 
subject = 'sub-04';

% choose dataset for headshape extraction
megfile = fullfile(subject,'meg',[subject,'_task-clicks_meg.fif']); 
%--------------------------------------------------------------------------

% make folder for data
%---------------------
new_dir = fullfile(settings.path2bids,'derivatives',subject,'layout');
if ~exist(new_dir, 'dir')
   mkdir(new_dir)
end


%% Prepare layout
%---------------
datapath = fullfile(settings.path2bids,megfile);
elec     = ft_convert_units(ft_read_sens(datapath,'senstype','eeg'),'mm'); 

cfg         = [];
cfg.elec    = elec;
layout_orig = ft_prepare_layout(cfg);

% Visualize layout
%-----------------
cfg        = [];
cfg.layout = layout_orig;
ft_layoutplot(cfg) 

% Save layout
%------------
save(fullfile(new_dir,[subject,'_layout-eeg_orig.mat']),'layout_orig') % in mm     

%% New labels
%-----------
label = {'Fp1';'Fp2';'F7';'F3';'FZ';'F4';'F8'; ...
         'T7';'C3';'Cz';'C4';'T8';'P7';'P3';'Pz';'P4';'P8'; ...
         'O1';'O2';'Fpz';'FC5';'FC1';'FC2';'FC6'; ...
         'CP5';'CP1';'CP2';'CP6';'PO5';'P0z';'PO6';'Oz'};    

elec.label        = label;
cfg               = [];
cfg.elec          = elec;
layout_relabelled = ft_prepare_layout(cfg);

% Save layout
%------------
save(fullfile(new_dir,[subject,'_layout-eeg_relabelled.mat']),'layout_relabelled') % in mm     

% Visualize layout
%-----------------
cfg        = [];
cfg.layout = layout_relabelled;
ft_layoutplot(cfg) 


