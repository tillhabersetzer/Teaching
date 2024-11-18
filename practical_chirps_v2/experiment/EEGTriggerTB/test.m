
ch = [12];
lch = length(ch);

ok = soundmexpro('init', ...     % command name
    'force',      1, ...
    'driver',     'ASIO Fireface', ...
    'samplerate', 48000, ...
    'output',     ch, ...
    'input',      [], ...
    'numbufs',    1, ...
    'autocleardata', 1, ...
    'ramplen', 0);


ok = soundmexpro('loadmem', ...      % command name
    'data', zeros(48,1), ...          % data vector
    'loopcount', 1 ...
    );

ok = soundmexpro('start','length',0);
if ~ok
    error(['error calling ''start''' ]);
end

for t = 1:1500
    t
    ok = soundmexpro('loadmem', ...      % command name
        'data',(2.^17-1)./(2.^24)  .*ones(24000,lch), ...          % data vector
        'loopcount', 1 ...
        );
    if ~ok
        error(['error calling ''loadmem''' ]);
    end
    pause(1)
end