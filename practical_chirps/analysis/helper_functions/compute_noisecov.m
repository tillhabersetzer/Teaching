function [avg_noise] = compute_noisecov(datapath,channel)
%--------------------------------------------------------------------------
% Till Habersetzer (28.01.2022)
% This script computes covariances matrices. It was initially designed to
% compute the noise covariance matrix given multiple empty room recordings.
% First, the recordings are appended and the covariance is computed for 
% each recording (e.g. trial) and the covariances are averaged afterwards 
% (see ft_ft_timelockanalysis).
%
% Input:
%   datapath: cell array with full path to recordings {1,F}
%   channel:  channeltype for computation 'meg','megmag','megplanar', etc.
%
% Output:
%   avg_noise: average of all recordings. The covariance is stored under
%              avg_noise.cov
%--------------------------------------------------------------------------

F     = length(datapath); % number of recordings, files
noise = cell(1,F);

% Load data
%----------
for f=1:F
    cfg              = [];
    cfg.dataset      = datapath{f};
    cfg.channel      = channel;
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0; % ensure that sensors are expressed in SI units
    noise{f}         = ft_preprocessing(cfg);
end

% Append data
%------------
if F>1
    cfg                = [];
    cfg.keepsampleinfo = 'no';
    noise              = ft_appenddata(cfg,noise{:});
else
    noise = noise{1};
end

% Compute covariance
%-------------------
cfg                  = [];
cfg.removemean       = 'yes'; % default for covariance computation
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
avg_noise            = ft_timelockanalysis(cfg,noise);

end
