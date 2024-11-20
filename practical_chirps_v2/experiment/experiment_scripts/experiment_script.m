function experiment_script(f,settings,name,stimulus,earphone)
%--------------------------------------------------------------------------
% parameters: 
%
% f:        figure handle
% settings: experimental settings 
% name:     subject code
% stimulus: 'click', 'up', 'down'
% earphone: earphones used
%--------------------------------------------------------------------------

%% Load Stimulus
%--------------------------------------------------------------------------
p = fileparts(mfilename('fullpath'));

% Stimulus Settings
switch earphone

    case {'tip300','sensimetrics'}
        audiofile = [stimulus,'.wav'];
        audio_dir = fullfile('audio','tip300');

        % Import audiofile
        [audio, fs] = audioread(fullfile(audio_dir,audiofile));
        % Check if mono
        if size(audio,2)~=1
            error('Mono file expected')
        end
        audio = [audio,audio];

    case 'sensimetrics_eq'
        audiofile = [stimulus,'_eq.wav'];
        audio_dir = fullfile('audio','sensimetrics');

        % Import audiofile
        [audio, fs] = audioread(fullfile(audio_dir,audiofile));
        % Check if stereo
        if size(audio,2)~=2
            error('Stereo file expected')
        end

end

% Check sampling frequency
if settings.samp_freq~=fs
    error("Sampling frequencies don't match (%i~=%i).",settings.samp_freq,fs)
end

%% Apply Calibration
%--------------------------------------------------------------------------
% calibration for transient based on peak-value
level_stim = 20*log10(max(abs(audio),[],1));

cal_val      = settings.calibration.(earphone).cal_val.(stimulus);
threshold    = settings.threshold.(earphone).(stimulus); % in dB (p-p) peSPL
target_level = threshold + settings.target_level_dB_SL; % in dB (p-p) peSPL

% Compute gain to achieve target_level
% level_stim + cal_val + gain_dB = target_level
% -> scaling = 10^(gain_dB/20)
gain_dB = target_level - level_stim - cal_val;
gain    = 10.^(gain_dB/20);  

cal_audio      = zeros(size(audio));
cal_audio(:,1) = audio(:,1) * gain(1); % left
cal_audio(:,2) = audio(:,2) * gain(2); % right
clear audio

%% Generate Trigger
%--------------------------------------------------------------------------
    
n_stim  = settings.n_stimuli;
jitter  = settings.jitter;
trig_id = settings.trig_id.(stimulus);     

switch stimulus
    case 'click'
        trig_len      = round(fs*0.01); % 10 ms
        trig_info     = EEGTrigID2info([trig_id],[16],'MEG'); % [8 8] distribution of 16 bits
        trig_amp      = trig_info.TrigWord;
        audio_trigger = trig_amp*ones(trig_len,1);

        rng('shuffle')
        % remove 10 ms because of extended stimulus
        jitter_duration_sec = (jitter(1)-0.01) + jitter(2)*rand(1,n_stim) ; % [0.34,0.39]s
        jitter_duration     = round(jitter_duration_sec*fs); % samples

        % Extend Stimulus to 10 ms
        %-------------------------
        cal_audio = vertcat(cal_audio,zeros(length(audio_trigger)-length(cal_audio),2));

    case {'up','down'}
        trig_len      = length(cal_audio);
        trig_info     = EEGTrigID2info([trig_id],16,'MEG'); % [8 8] distribution of 16 bits
        trig_amp      = trig_info.TrigWord;
        audio_trigger = trig_amp*ones(trig_len,1);

        rng('shuffle')
        jitter_duration_sec = jitter(1) + jitter(2)*rand(1,n_stim) ; % [0.35,0.40]s
        jitter_duration     = round(jitter_duration_sec*fs); % samples
end

%% Stimulus Presentation with SoundMexPro 
%--------------------------------------------------------------------------
asio_driver = settings.asio_driver;
if strcmp(asio_driver,'ASIO Fireface USB')
    output_channel = [0 1 8];
elseif strcmp(asio_driver,'Focusrite USB ASIO')
    output_channel = [0 1 2];
else 
    error('Unexpected Asio Driver: %s!',asio_driver)
end

ok = soundmexpro('exit');
if ~ok, error('error calling ''exit'' '); end             

ok = soundmexpro('init', ...       % command name
    'samplerate',   settings.samp_freq, ...
    'force',        1, ...         % exit called internally before init
    'driver',       asio_driver, ...
    'output',       output_channel, ...   % soundcard channels [left,right,trigger]
    'input',        -1, ...        % no input channel used
    'track',         3, ...        % [stim stim trigger] [left right trigger]
    'autocleardata', 1);           % audio data already played completely should be cleared from memory automatically on next data loading command
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('trackmap','track', [0 1 2]); % virtual channel (first track -> channel 0, second track -> channel 1, third track -> channel 3)        
if ~ok, error('error calling ''trackmap'' '); end             
 
ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

% Play Zeros
%--------------------------------------------------------------------------

ok = soundmexpro('loadmem', ...
     'data',zeros(1*fs,3), ...
     'track',[0 1 2], ...
     'loopcount', 1);
if ~ok,  error('error calling ''loadmem''' ); end

% Play Stimuli
%--------------------------------------------------------------------------

% progress bar
d = uiprogressdlg(f,...
                  'Title','Tracks loaded',...
                  'Icon','info',...
                  'ShowPercentage','on',...
                  'Cancelable','on');

% Length
% L = sum(JitterDuration) + Nstim*length(Stim) + SampFreq*2;

% Counter
tracks_loaded = 0;
flipsign      = zeros(1,n_stim);

% Preload data
fsign = 1;
for n = 1:10
    ok = soundmexpro('loadmem', ...
     'data',vertcat([fsign*cal_audio, audio_trigger],zeros(jitter_duration(tracks_loaded+1),3)), ...
     'track',[0 1 2], ...
     'loopcount', 1);
    if ~ok,  error('error calling ''loadmem''' ); end

    tracks_loaded           = tracks_loaded + 1;
    d.Value                 = tracks_loaded/n_stim;
    d.Message               = sprintf('Track %i of %i loaded.',tracks_loaded, n_stim); 
    flipsign(tracks_loaded) = fsign;
    fsign                   = -fsign;
    pause(0.01)
end

ok = soundmexpro('start','length', 0);             
if ~ok,  error('error calling ''start''' ); end

while tracks_loaded<n_stim

    [ok, trackload] = soundmexpro('trackload');
    if ok ~= 1
        error('error calling ''trackload''');
    end

    if d.CancelRequested
         break
    end
    
    if (trackload(1) < 10)
        
        ok = soundmexpro('loadmem', ...            
            'data', vertcat([fsign*cal_audio, audio_trigger],zeros(jitter_duration(tracks_loaded+1),3)), ...     
            'track',[0 1 2], ...                   
            'loopcount', 1);                       
        if ~ok  
            error('error calling ''loadmem'' ');
        end
        tracks_loaded           = tracks_loaded + 1;
        d.Value                 = tracks_loaded/n_stim;
        d.Message               = sprintf('Track %i of %i loaded.',tracks_loaded, n_stim); 
        flipsign(tracks_loaded) = fsign;
        fsign                   = -fsign;
        pause(0.01)
    end  

end

pause(5) % pause is needed for preloaded buffer: 10 files

ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end   

disp('Playback of experimental_script.m finished successfully.')

% Save used Settings
%--------------------------------------------------------------------------
answer = questdlg('Should the measurement data be saved?','Saving','Yes','No','Yes');
switch answer
    case 'Yes'
        pth = fullfile(p,'..','results',name);
        if ~exist(pth,'dir')
            mkdir(pth)
        end
        results.settings            = settings;
        results.measdate            = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % letztes Update der Einstellungen
        results.jitter              = jitter_duration_sec+0.01;
        results.flipsign            = flipsign;
        results.used_cal_val        = cal_val;
        results.used_threshold      = threshold;
        results.target_level        = target_level;
        results.level_stim          = level_stim;
        results.selection.subjectid = name;
        results.selection.earphone  = earphone;
        results.selection.stimulus  = stimulus;
        results.selection.earphone  = earphone;
        switch stimulus
            case 'click'
                results.jitter = jitter_duration_sec+0.01;
            case {'up','down'}
                results.jitter = jitter_duration_sec;
        end
          
        % Avoid overwriting of existing file
        in_file  = sprintf('%s_stim-%s_%s_results.mat',name,stimulus,earphone); 
        in_path  = pth;
        out_file = avoidOverwrite(in_file,in_path,2,0);
        save(fullfile(pth,out_file),'results')
    case 'No'
end

end