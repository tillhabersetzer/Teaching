close all
clear all
clc

time = 0:1/44100:1;
sig  = 0.01*sin(2*pi*1000*time);

ASIODriver = 'Focusrite USB ASIO';

ok = soundmexpro('exit');
if ~ok, error('error calling ''exit'' '); end             

ok = soundmexpro('init', ...       % command name
    'samplerate',   44100, ...
    'force',        1, ...         % exit called internally before init
    'driver',       ASIODriver, ...
    'output',       [0 1], ...   % soundcard channels [left,right,trigger]
    'input',        -1, ...        % no input channel used
    'track',         2, ...        % [stim stim trigger] [left right trigger]
    'autocleardata', 1);           % audio data already played completely should be cleared from memory automatically on next data loading command
if ~ok, error('cannot initialize soundmexpro, error calling ''init'' '); end

% Preload data
fsign = 1;
for n = 1:10
    ok = soundmexpro('loadmem', ...
     'data',[sig',sig'], ...
     'track',[0 1], ...
     'loopcount', 1);
    if ~ok,  error('error calling ''loadmem''' ); end
end

L = 10*length(sig);

ok = soundmexpro('start','length', L);             
if ~ok,  error('error calling ''start''' ); end

% ok = soundmexpro('stop'); if ~ok,  error('error calling ''stop'' '); end    
% ok = soundmexpro('exit'); if ~ok,  error('error calling ''exit'' '); end   