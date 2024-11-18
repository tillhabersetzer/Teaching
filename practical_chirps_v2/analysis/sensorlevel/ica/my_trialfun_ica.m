function [trl, event] = my_trialfun_ica(cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trialfunction for ica 
% Selects the interval between the first and the last appearance of the
% trigger with specified eventtype. The interval is than 
% segmented into segments with specified lenghts in seconds. Thereby the
% seconds are defined by the sampling frequency.
%
% cfg.trialdef.eventtype = 'STI101'; chirps + story interval - highly
% recommmend because it's less vulnerable to trigger errors
% cfg.trialdef.eventtype = 'STI002'; only chirp interval
% cfg.trialdef.eventtype = 'STI001'; only story interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

smp = [event.sample];
typ = {event.type};

% fullfill eventtype 
idx = find(strcmp(typ,cfg.trialdef.eventtype));

% interval length in samples
stepsize = round(cfg.trialdef.triallength * hdr.Fs); % in samples
first    = smp(idx(1));
last     = smp(idx(end));
trlbegin = (first:stepsize:last)';

if trlbegin(end) == last
    trlbegin(end) = [];  
end

trlend = [trlbegin(2:end)-1;last];
trl    = [trlbegin trlend trlbegin-first];

end

