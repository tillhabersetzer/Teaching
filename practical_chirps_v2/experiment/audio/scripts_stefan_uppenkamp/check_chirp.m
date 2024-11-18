%--------------------------------------------------------------------------
% Till Habersetzer, 22.10.2024
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%--------------------------------------------------------------------------

close all
clearvars
clc

%% Compute ffT
%--------------------------------------------------------------------------

% optionally 
% cancatenate 100 Stimuli and Compute spectrum
% concatenate_stim = true;
concatenate_stim = false;

% Load audio signals
%-------------------
% audionames = {'click','up','down','up_dau'};
audionames = {'click','up','down'};
S          = length(audionames);
sig        = cell(1,S);
audio      = cell(1,S);
fsamp      = cell(1,S);
time       = cell(1,S);  
freqs      = cell(1,S);
fft_spec   = cell(1,S);
n_fft      = zeros(1,S);
n_stim     = 100;
gap        = 0.0;
polarity   = 1;

figure
for s = 1:S
    [audio{s},fsamp{s}] = audioread(fullfile('final_stimuli',[audionames{s},'.wav']));  
    % [audio{s},fsamp{s}] = audioread(fullfile('..','tip300',[audionames{s},'.wav']));  
    time{s}             = (0:length(audio{s})-1)*1000/fsamp{s}; % in ms

    if concatenate_stim
        sig{s} = [];
        for nidx = 1:n_stim
            switch audionames{s}
                case 'click'
                    % add 10 ms zeros anyway
                    sig{s} = [sig{s};audio{s};zeros(fsamp{s}*0.01,1);zeros(fsamp{s}*gap,1)];
                case {'up','down'}
                    if mod(nidx,2) 
                        sig{s} = [sig{s};polarity*audio{s};zeros(fsamp{s}*gap,1)];
                    else % flip even signals
                        % sig{s} = [sig{s};polarity*flipud(audio{s});zeros(fsamp{s}*gap,1)];
                        sig{s} = [sig{s};polarity*audio{s};zeros(fsamp{s}*gap,1)];
                    end
            end

            % sig{s} = [sig{s};audio{s}] ;
            % polarity = -polarity;
        end
        n_fft(s) = 2^nextpow2(length(sig{s}));

    else
        switch audionames{s}
            case 'click'
                % add 10 ms zeros anyway
                sig{s} = [audio{s};zeros(fsamp{s}*0.01,1)];
            otherwise
                sig{s} = audio{s};
        end
        n_fft(s) = 2^nextpow2(n_stim*(length(sig{s})+fsamp{s}*gap));
    end

    subplot(S,1,s)
    plot(sig{s})
    title(audionames{s})

    % Compute spectrum
    %-----------------
    % L           = length(sig{s});
    L           = n_fft(s);
    Y           = fft(sig{s},n_fft(s));
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

axs = cell(1,S);
figure
for s = 1:S
    subplot(S,2,2*s-1); % odd numbers
    plot(time{s},audio{s},'-x')
    xlabel('t / ms');
    title(audionames{s})
    grid('on')
    
    axs{s} = subplot(S,2,2*s);
    plot(freqs{s}/1000,db(fft_spec{s},'power')) % dB: 10*log10(x)
    hold on
    xlabel("f / kHz")
    ylabel("|P1(f)|^2 / dB")
    xlim([0.001,15])
    % xlim([f1*0.001/10,2*f2*0.001])
    title(audionames{s})
    grid('on')
end
set([axs{:}],'xtick', [10,100,1000,5000,10000,15000]/1000, 'xticklabel',{'0.01' '0.1' '1' '5' '10' '15'})
set([axs{:}],'XScale', 'log');

% s=2;
% soundsc(sig{s},fsamp{s})
