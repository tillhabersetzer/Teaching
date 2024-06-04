function noise = give_noise(filenames,subjectdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenates multiple noise files
% filenames: cell array with noise measurements filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N          = length(filenames);
noise      = cell(1,N);

for n = 1:N
    cfg                = [];
    cfg.dataset        = fullfile(subjectdata.rawdatadir,[filenames{n},'.fif']);
    cfg.channel        = 'meg'; 
    cfg.continuous     = 'yes';
    cfg.coilaccuracy   = 0;            % ensure that sensors are expressed in SI units
    cfg.demean         = 'yes';        % important for covariance estimate
    cfg.baselinewindow = 'all';
    noise{n}           = ft_preprocessing(cfg);   
end

hdr                = noise{1}.hdr;
cfg                = [];
cfg.keepsampleinfo = 'no';
noise              = ft_appenddata(cfg,noise{:});
noise.hdr          = hdr;

noisi = [];
timi  = [];
for t = 1:length(noise.trial)
    noisi    = [noisi,noise.trial{t}];
    timi_add = noise.time{t};
    if t>1  
        timi_add = (timi(end)+1/noise.fsample)+timi_add; % add time offset
    end
    timi = [timi,timi_add];
end
noise.trial = {noisi};
noise.time  = {timi};

end