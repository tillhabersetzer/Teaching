close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
addpath(fullfile('..','..'))
eval('main_settings')

addpath(fullfile('..','..','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
subjects = [1,2];

% Choose timewindow for plotting
time2plot = [-50 300];
%--------------------------------------------------------------------------

%% Load data
%--------------------------------------------------------------------------
% Arrange data: Sx3: subjects x conditions
S            = length(subjects);
sources      = [];
sources_mean = [];
sources_svd  = [];
dippos_nosym = cell(1,S);

counter = 1;
for subidx=subjects
    
    subject = ['sub-',num2str(subidx,'%02d')];
    data    = importdata(fullfile(settings.path2project,'derivatives',subject,'sourcelevel',[subject,'_dipolefits-N19mP30m.mat']));
   
    sources      = vertcat(sources,{data.source_vec{1}.dip.mom,data.source_vec{2}.dip.mom,data.source_vec{3}.dip.mom});
    sources_mean = vertcat(sources_mean,data.source_sca_mean);
    sources_svd  = vertcat(sources_svd,data.source_sca_svd);

    dippos_nosym{counter} = data.source_pooled_nosym.dip.pos;
    counter = counter +1;
end

time = data.source_vec{1}.time*1000;
clear data

% Rescale units - optional for plotting (Am -> nAm)
%--------------------------------------------------
for subidx=1:S
    for cidx = 1:3
        sources{subidx,cidx}      = 10^9*sources{subidx,cidx};
        sources_mean{subidx,cidx} = 10^9*sources_mean{subidx,cidx};
        sources_svd{subidx,cidx}  = 10^9*sources_svd{subidx,cidx};
    end
end

%% Compute Grand averages
%--------------------------------------------------------------------------

gavg_sources      = cell(1,3);
gavg_sources_mean = cell(1,3);
gavg_sources_svd  = cell(1,3);

for cidx=1:3 % loop over conditions
    for subidx=1:S % loop over subjects

        % Mapping between dipole locations and hemispheres
        %-------------------------------------------------
        mapping = check_diploc(dippos_nosym{subidx});
        % order right and left hemispere, 1st row: left, 2nd row: right
        sources_mean{subidx,cidx} = sources_mean{subidx,cidx}(mapping,:);
        sources_svd{subidx,cidx}  = sources_svd{subidx,cidx}(mapping,:);

        timecourse_mean(:,:,subidx) = sources_mean{subidx,cidx};
        timecourse_svd(:,:,subidx)  = sources_svd{subidx,cidx};

        % for xyz dipole moments
        if mapping(1)==1
            sources_mean{subidx,cidx} = sources{subidx,cidx}([1:3,4:6],:);
        else
            sources_mean{subidx,cidx} = sources{subidx,cidx}([4:6,1:3],:);            
        end
        timecourse(:,:,subidx) = sources{subidx,cidx};

    end
    gavg_sources{cidx}      = mean(timecourse,3); % mean per condition
    gavg_sources_mean{cidx} = mean(timecourse_mean,3); % mean per condition
    gavg_sources_svd{cidx}  = mean(timecourse_svd,3); % mean per condition
end
clear timecourse_mean timecourse_svd timecourse

%% All moments, xyz-directions and both hemispheres
%--------------------------------------------------------------------------

name = {'left hemisphere','right hemisphere'};
% latency correction

idxs      = [];
idxs(1,:) = 1:3; % left
idxs(2,:) = 4:6; % right

figidxs = [1,3,5;2,4,6];

figure
for hidx=1:2
    
    idx    = idxs(hidx,:);
    figidx = figidxs(hidx,:);

    mini = min([gavg_sources{1}(idx,:),gavg_sources{2}(idx,:),gavg_sources{3}(idx,:)],[],'all');
    maxi = max([gavg_sources{1}(idx,:),gavg_sources{2}(idx,:),gavg_sources{3}(idx,:)],[],'all');
  
    axisvec = horzcat(time2plot,[mini,maxi]);
    
    subplot(3,2,figidx(1)); 
    plot(time, gavg_sources{1}(idx,:), '-')
    if hidx==1; ylabel('dipole moment / nAm'); end
    legend({'x', 'y', 'z'});
    axis(axisvec) 
    grid on
    title(['click ',name{hidx}])
    
    subplot(3,2,figidx(2)); 
    plot(time, gavg_sources{2}(idx,:), '-')
    if hidx==1; ylabel('dipole moment / nAm'); end
    legend({'x', 'y', 'z'});
    axis(axisvec)
    grid on
    title(['upchirp ',name{hidx}])
    
    subplot(3,2,figidx(3)); 
    plot(time, gavg_sources{3}(idx,:), '-')
    xlabel('t / ms')
    if hidx==1; ylabel('dipole moment / nAm'); end
    legend({'x', 'y', 'z'});
    axis(axisvec)
    grid on
    title(['downchirp ',name{hidx}])
end

%% Visualize fixed dipole (fixed orientation)
%--------------------------------------------------------------------------
% mean dipolmoment orientation has been used as orientation constraint 

name = {'left hemisphere','right hemisphere'};
% latency correction

figure
for hidx=1:2

    mini = min([gavg_sources_mean{1}(hidx,:),gavg_sources_mean{2}(hidx,:),gavg_sources_mean{3}(hidx,:)],[],'all');
    maxi = max([gavg_sources_mean{1}(hidx,:),gavg_sources_mean{2}(hidx,:),gavg_sources_mean{3}(hidx,:)],[],'all');

    axisvec = horzcat(time2plot,[mini,maxi]);
    
    subplot(2,1,hidx); 
    hold on
    plot(time, gavg_sources_mean{1}(hidx,:), '-')
    plot(time, gavg_sources_mean{2}(hidx,:), '-')
    plot(time, gavg_sources_mean{3}(hidx,:), '-')
    legend({'click', 'upchirp', 'downchirp'});
    axis(axisvec)
    grid on
    title(name{hidx})
    if hidx==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm')
end
sgtitle('gavg: fixed oriented dipole via mean orientation')

%% Visualize fixed dipole (fixed orientation)
%--------------------------------------------------------------------------
% maximum variance dipolmoment orientation via SVDhas been used as 
% orientation constraint

name = {'left hemisphere','right hemisphere'};
% latency correction

figure
for hidx=1:2

    mini = min([gavg_sources_svd{1}(hidx,:),gavg_sources_svd{2}(hidx,:),gavg_sources_svd{3}(hidx,:)],[],'all');
    maxi = max([gavg_sources_svd{1}(hidx,:),gavg_sources_svd{2}(hidx,:),gavg_sources_svd{3}(hidx,:)],[],'all');

    axisvec = horzcat(time2plot,[mini,maxi]);
    
    subplot(2,1,hidx); 
    hold on
    plot(time, gavg_sources_svd{1}(hidx,:), '-')
    plot(time, gavg_sources_svd{2}(hidx,:), '-')
    plot(time, gavg_sources_svd{3}(hidx,:), '-')
    legend({'click', 'upchirp', 'downchirp'});
    axis(axisvec)
    grid on
    title(name{hidx})
    if hidx==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm')
end
sgtitle('gavg: fixed oriented dipole via svd orientation')
