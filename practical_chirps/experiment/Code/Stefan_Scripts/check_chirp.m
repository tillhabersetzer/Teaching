%--------------------------------------------------------------------------
% Till Habersetzer, 22.10.2024
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%--------------------------------------------------------------------------

close all
clear all
clc

%% Compute ffT
%--------------------------------------------------------------------------

% Load audio signals
%-------------------
audionames = {'click','up','down'};
S          = length(audionames);
sig        = cell(1,S);
fsamp      = cell(1,S);
time       = cell(1,S);  
freqs      = cell(1,S);
fft_spec   = cell(1,S);

n_fft = 4096;

for s = 1:S
    [sig{s},fsamp{s}] = audioread(fullfile('final_stimuli',[audionames{s},'.wav']));  
    time{s}           = (0:length(sig{s})-1)*1000/fsamp{s}; % in ms

    % Compute spectrum
    %-----------------
    % L           = length(sig{s});
    L           = n_fft;
    Y           = fft(sig{s},n_fft);
    P           = abs(Y/L).^2; % normalize and power
    P           = P(1:L/2+1);
    freqs{s}    = fsamp{s}*(0:(L/2))/L;
    fft_spec{s} = P;
    clear L Y P
end

if isequal(fsamp{:})
    fs = fsamp{1};
else
    error('fs not consistent')
end

figure
for s = 1:S
    subplot(S,2,2*s-1); % odd numbers
    plot(time{s},sig{s},'-x')
    xlabel('t / ms');
    title(audionames{s})
    
    subplot(S,2,2*s)
    semilogx(freqs{s}/1000,db(fft_spec{s},'power')) % dB: 10*log10(x)
    hold on
    xlabel("f / kHz")
    ylabel("|P1(f)|^2 / dB")
    %xlim([f1*0.001/10,2*f2*0.001])
    title(audionames{s})
end

%% Old checkup
%--------------------------------------------------------------------------

% check_wav   = 1;
% create_wav  = 0;
% concat_stim = 1;
% 
% if create_wav
% %     fstart      = 100;
% %     fstop       = 10000;
% %     fsamp       = 44100;
% %     period      = 0.01;
% %     a           =  2;
% %     duration    = 0.5;
% %     scalefactor = 30000;
% %     make_chirps(fstart, fstop, fsamp, period, a, duration, scalefactor);
% %     make_chirps_flat(fstart, fstop, fsamp, period, a, duration, scalefactor);
% %     fs     = fsamp;
%     f1     = 100;
%     f2     = 10000;
%     fs     = 44100;
%     len    = 0;
%     alpha  = 3;
%     sig{1} = bmchirp(f1,f2,fs,len);
%     sig{2} = alchirp(f1,f2,fs,len,alpha);
%     S      = length(sig);
%     audionames = {'bmchirp','alchirp'};
% end
% 
% if check_wav
%     audionames = {'click','gdn22','gdn30','gup22','gup30','sound1','up','down'};
%     S          = length(audionames);
%     sig        = cell(1,S);
%     fsamp      = cell(1,S);
% 
%     for s = 1:S
%         [sig{s},fsamp{s}] = audioread([audionames{s} ,'.wav']);   
%     end
% 
%     if isequal(fsamp{:})
%         fs = fsamp{1};
%     else
%         error('fs not consistent')
%     end
% 
% end
% 
% if concat_stim
% 
%     N = 10;
% 
%     for s = 1:S
%         L      = length(sig{s});
%         sig{s} = sig{s}.*hann(L)';
%         sig{s} = repmat(sig{s},1,N);
%     end
% 
% end
% 
% pxx_pw  = cell(1,S); % pwelch
% freq_pw = cell(1,S);
% pxx_pd  = cell(1,S); % periodogram
% freq_pd = cell(1,S);
% time    = cell(1,S);
% 
% for s = 1:S
%     [pxx_pw{s},freq_pw{s}] = pwelch(sig{s},[],[],[],fs);
%     [pxx_pd{s},freq_pd{s}] = periodogram(sig{s},[],[],fs);    
%     time{s}                = (0:length(sig{s})-1)*1000/fs; % in ms
% end
% 
% figure
% for s = 1:S
% subplot(S,2,2*s-1)
% plot(time{s},sig{s},'-x')
% xlabel('t / ms');
% title(audionames{s})
% 
% subplot(S,2,2*s)
% plot(freq_pw{s}/1000,db(pxx_pw{s},'power'))
% hold on
% plot(freq_pd{s}/1000,db(pxx_pd{s},'power'))
% xlabel('f / kHz')
% %xlim([f1*0.001/10,2*f2*0.001])
% title(audionames{s})
% legend({'pwelch','periodogram'})
% end
% 




