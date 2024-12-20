function experimental_script(f,Settings,name,type)
%--------------------------------------------------------------------------
% parameters: 
%
% f:        figure handle
% Settings: experimental settings 
% name:     subject code
% type:     'click', 'up', 'down'
%--------------------------------------------------------------------------

% Load Stimulus
%==========================================================================
switch type
    case 'click'
        filename = 'click_rmseq.wav';
        TrigID   = Settings.TrigID.click;            
    case 'up'
        filename = 'up_rmseq.wav';
        TrigID   = Settings.TrigID.upchirp;
    case 'down'
        filename = 'down_rmseq.wav';
        TrigID   = Settings.TrigID.downchirp;      
end

p               = fileparts(mfilename('fullpath'));
[stim,SampFreq] = audioread(fullfile(p,'..','Stimuli','rms_equalized',filename));

if Settings.SampFreq~=SampFreq
    error("Sampling frequencies don't match (%i~=%i).",Settings.SampFreq,SampFreq)
end

% Apply Calibration
%==========================================================================
% AttdB + RMSdBStim + CalVal = TargetLevel
% scaling = 10^(AttdB/20)

AttFact_stim = 10.^(Settings.Calibration.AttdB/20);
Stim         = zeros(length(stim),2);
Stim(:,1)    = stim * AttFact_stim(1); % left
Stim(:,2)    = stim * AttFact_stim(2); % right

% Generate Trigger
%==========================================================================
    
Nstim  = Settings.Nstimuli;
Jitter = Settings.Jitter;

switch type
    case 'click'
        TrigLen      = round(SampFreq*0.01); % 10 ms
        TrigInfo     = EEGTrigID2info(TrigID,16,'MEG'); % [8 8] distribution of 16 bits
        TrigAmp      = TrigInfo.TrigWord;
        stim_trigger = TrigAmp*ones(TrigLen,1);

        rng('shuffle')
        % remove 10 ms because of extended stimulus
        JitterDuration_sec = (Jitter(1)-0.01) + Jitter(2)*rand(1,Nstim) ; % [0.34,0.39]s
        JitterDuration     = round(JitterDuration_sec*SampFreq); % samples

        % Extend Stimulus to 10 ms
        %=========================
        Stim = vertcat(Stim,zeros(length(stim_trigger)-length(Stim),2));

    case {'up','down'}
        TrigLen      = length(Stim);
        TrigInfo     = EEGTrigID2info(TrigID,16,'MEG'); % [8 8] distribution of 16 bits
        TrigAmp      = TrigInfo.TrigWord;
        stim_trigger = TrigAmp*ones(TrigLen,1);

        rng('shuffle')
        JitterDuration_sec = Jitter(1) + Jitter(2)*rand(1,Nstim) ; % [0.35,0.40]s
        JitterDuration     = round(JitterDuration_sec*SampFreq); % samples
end

% Stimulus Presentation with SoundMexPro 
%==========================================================================
ASIODriver = Settings.Driver;

ok = soundmexpro('exit');
if ~ok, error('error calling ''exit'' '); end             

ok = soundmexpro('init', ...       % command name
    'samplerate',   SampFreq, ...
    'force',        1, ...         % exit called internally before init
    'driver',       ASIODriver, ...
    'output',       [0 1 8], ...   % soundcard channels [left,right,trigger]
    'input',        -1, ...        % no input channel used
    'track',         3, ...        % [stim stim trigger] [left right trigger]
    'autocleardata', 1);           % audio data already played completely should be cleared from memory automatically on next data loading command
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

ok = soundmexpro('trackmap','track', [0 1 2]); % virtual channel (first track -> channel 0, second track -> channel 1, third track -> channel 3)        
if ~ok, error('error calling ''trackmap'' '); end             
 
ok = soundmexpro('cleardata');
if ~ok, error('error calling ''cleardata'' '); end

% Play Zeros
%==========================================================================

% switch on debug saving
% if 1 ~= soundmexpro('debugsave', ...    % command name
%         'value', [1], ...               % enable/disable flag (here: enable)
%         'channel', [2] ...              % channel(s) where to apply value
%         );
%     error(['error calling ''debugsave''' error_loc(dbstack)]);
% end

ok = soundmexpro('loadmem', ...
     'data',zeros(1*SampFreq,3), ...
     'track',[0 1 2], ...
     'loopcount', 1);
if ~ok,  error('error calling ''loadmem''' ); end

% Play Stimuli
%==========================================================================

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
flipsign      = zeros(1,Nstim);

% Preload data
fsign = 1;
for n = 1:10
    ok = soundmexpro('loadmem', ...
     'data',vertcat([fsign*Stim, stim_trigger],zeros(JitterDuration(tracks_loaded+1),3)), ...
     'track',[0 1 2], ...
     'loopcount', 1);
    if ~ok,  error('error calling ''loadmem''' ); end

    tracks_loaded           = tracks_loaded + 1;
    d.Value                 = tracks_loaded/Nstim;
    d.Message               = ['Track ',num2str(tracks_loaded),' of ',num2str(Nstim),' loaded.']; 
    flipsign(tracks_loaded) = fsign;
    fsign                   = -fsign;
    pause(0.01)
end

ok = soundmexpro('start','length', 0);             
if ~ok,  error('error calling ''start''' ); end

while tracks_loaded<Nstim

    [ok, trackload] = soundmexpro('trackload');
    if ok ~= 1
        error('error calling ''trackload''');
    end

    if d.CancelRequested
         break
    end
    
    if (trackload(1) < 10)
        
        ok = soundmexpro('loadmem', ...            
            'data', vertcat([fsign*Stim, stim_trigger],zeros(JitterDuration(tracks_loaded+1),3)), ...     
            'track',[0 1 2], ...                   
            'loopcount', 1);                       
        if ~ok  
            error('error calling ''loadmem'' ');
        end
        tracks_loaded           = tracks_loaded + 1;
        d.Value                 = tracks_loaded/Nstim;
        d.Message               = ['Track ',num2str(tracks_loaded),' of ',num2str(Nstim),' loaded.'];  
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
        pth = fullfile(p,'..','..','Results',name);
        if ~exist(pth,'dir')
            mkdir(pth)
        end
        info.Settings = Settings;
        info.Measdate = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')); % letztes Update der Einstellungen
        info.flipsign = flipsign;
        switch type
            case 'click'
                info.Jitter = JitterDuration_sec+0.01;
            case {'up','down'}
                info.Jitter = JitterDuration_sec;
        end
          
        % Avoid overwriting of existing file
        inFile  = [name,'_',type,'_info.mat'];
        inPath  = pth;
        outFile = avoidOverwrite(inFile,inPath,2,0);
        save(fullfile(pth,outFile),'info')
    case 'No'
end

end