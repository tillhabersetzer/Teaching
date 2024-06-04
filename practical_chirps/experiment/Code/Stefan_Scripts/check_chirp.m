close all; clear all; clc

check_wav  = 1;
create_wav = 0;

if create_wav
%     fstart      = 100;
%     fstop       = 10000;
%     fsamp       = 44100;
%     period      = 0.01;
%     a           =  2;
%     duration    = 0.5;
%     scalefactor = 30000;
%     make_chirps(fstart, fstop, fsamp, period, a, duration, scalefactor);
%     make_chirps_flat(fstart, fstop, fsamp, period, a, duration, scalefactor);
%     fs     = fsamp;
    f1     = 100;
    f2     = 10000;
    fs     = 44100;
    len    = 0;
    alpha  = 3;
    sig{1} = bmchirp(f1,f2,fs,len);
    sig{2} = alchirp(f1,f2,fs,len,alpha);
    S      = length(sig);
    audionames = {'bmchirp','alchirp'};
end

if check_wav
    audionames = {'click','gdn22','gdn30','gup22','gup30','sound1','up','down'};
    S          = length(audionames);
    sig        = cell(1,S);
    fsamp      = cell(1,S);

    for s = 1:S
        [sig{s},fsamp{s}] = audioread([audionames{s} ,'.wav']);   
    end
    
    if isequal(fsamp{:})
        fs = fsamp{1};
    else
        error('fs not consistent')
    end

end

pxx_pw  = cell(1,S); % pwelch
freq_pw = cell(1,S);
pxx_pd  = cell(1,S); % periodogram
freq_pd = cell(1,S);
time    = cell(1,S);

for s = 1:S
[pxx_pw{s},freq_pw{s}] = pwelch(sig{s},[],[],[],fs);
[pxx_pd{s},freq_pd{s}] = periodogram(sig{s},[],[],fs);    
time{s}                = (0:length(sig{s})-1)*1000/fs; % in ms
end

figure
for s = 1:S
subplot(S,2,2*s-1)
plot(time{s},sig{s},'-x')
xlabel('t / ms');
title(audionames{s})

subplot(S,2,2*s)
plot(freq_pw{s}/1000,db(pxx_pw{s},'power'))
hold on
plot(freq_pd{s}/1000,db(pxx_pd{s},'power'))
xlabel('f / kHz')
%xlim([f1*0.001/10,2*f2*0.001])
title(audionames{s})
legend({'pwelch','periodogram'})
end
