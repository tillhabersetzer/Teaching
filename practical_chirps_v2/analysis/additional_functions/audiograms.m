%--------------------------------------------------------------------------
% Till Habersetzer, 11.11.2024
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%--------------------------------------------------------------------------

close all 
clearvars 
clc

freq = [125 250 500 750 1000 1500 2000 3000 4000 6000 8000]; % Hz

% contains Hörpegel in dB
audiogram_r = [];
audiogram_l = [];

name        = {};

% participants
%-------------
% 1.) sub-01 / IN06ED16 - 05.11.24

name{1}          = 'sub-01';
audiogram_r(1,:) = [0,0,-5,-5,-5,0,0,0,5,10,20];
audiogram_l(1,:) = [-5,-5,-5,-5,-5,0,0,0,5,10,5];

% 2.) sub-02 / ST08LF13 - 05.11.24

name{2}          = 'sub-02';
audiogram_r(2,:) = [5,0,0,-5,0,0,0,0,0,5,0];
audiogram_l(2,:) = [10,5,5,10,5,5,0,-5,5,5,-5];

% 3.) sub-03 / KO07KL26 - 06.01.20

name{3}          = 'sub-03';
audiogram_r(3,:) = [5,5,5,5,0,-5,-5,-10,-10,-10,-10];
audiogram_l(3,:) = [5,5,0,5,5,5,5,-5,-10,-10,-10];

%--------------------------------------------------------------------------
Audiograms_L = array2table([freq; audiogram_l],'RowNames',['freq';name']);
Audiograms_R = array2table([freq; audiogram_r],'RowNames',['freq';name']);

%% Plot Audiograms
%--------------------------------------------------------------------------
N = length(name);
figure('Name','Audiogramme')

axs = cell(1,2);

axs{1} = subplot(1,2,1);
for i = 1:N
    plot(freq/1000,audiogram_r(i,:),'-o');
    hold on
end
legend(name,'Location','Southwest');
title('Right')
xlabel('f / lHz')
ylabel(' Level / dB')
grid on

axs{2} = subplot(1,2,2);
for i = 1:N
    plot(freq/1000,audiogram_l(i,:),'-o');
    hold on
end
legend(name,'Location','Southwest');
title('Left')
xlabel('f / kHz')
ylabel(' Level / dB')
grid on

set([axs{:}],'xtick', [125 250 500 1000 2000 4000 8000]/1000, 'xticklabel',{'0.125' '0.25' '0.5' '1' '2' '4', '8'})
set([axs{:}],'XScale', 'log');
set([axs{:}], 'YDir','reverse','ylim',[-10,30])

%% Convert data into "bids compatible" format
%--------------------------------------------------------------------------

rawdata_path = fullfile('..','..','rawdata');

% Define the list of subject IDs
for subidx = 1:3

    subject   = sprintf('sub-%02d', subidx);
    audiogram = array2table([freq', audiogram_l(subidx,:)', audiogram_r(subidx,:)'],'VariableNames',{'frequency','hearing_threshold_left','hearing_threshold_right'});
   
    % Save the table as a .tsv file
    filename = fullfile(rawdata_path,subject,'beh', [subject '_task-audiogram_beh.tsv']);
    writetable(audiogram, filename, 'FileType', 'text', 'Delimiter', '\t');
end
