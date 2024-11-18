function [scalar_dipmom] = constrained_dipolfitting(cfg,data)
%--------------------------------------------------------------------------
% Till Habersetzer (10.02.2022)
% This function performs dipolfitting of two dipols ins case of predefined 
% dipolmoments and orientations.
% The current fieldtrip version doesnt support dipolfitting with fixed
% orientations and consequently returns freely oriented dipols with xyz-
% directions for the dipolmoment. 
% This functions computs the scalar dipolmoments for 2 dipoles with given
% location and orientation. 
% The orientation is computed as a mean over a precalculated dipolmoment
% timecourse (3,T) or can be given in cfg.ori.
% The solution of the overconstrained linear equation can be found in:
% Brain signals: Physics and Mathematics of MEG and EEG, p.97, eq. (5.23)
%
%
% Output:
%   scalar_dipmom: scalar dipolmoment for 2 dipols (2,T)
%
% Input: 
%   cfg.dippos:    dipol positions (2,3), each row contains 1 dipole position
%   cfg.unit:      unit of dipol positions
%   cfg.headmodel: precomputed headmodel
%   cfg.channel:   channeltype for leadfield computation 'meg' 'megmag'
%                  'megplanar' 'eeg' etc.
%   cfg.noisecov:   optional, if provided: sphere the date using the noise 
%                   covariance matrix supplied
%   data:           averaged timecourse of data, usually computed via
%                   ft_timelockanalysis
%
%   Use either cfg.ori or cfg.dipmom
%   cfg.ori:       Assumed dipol orientations (3,2). Each column
%                  comprises one orientation.
%   cfg.dipmom     Only used if no orientation is provided. The
%                  orientations are calculated as a mean over a freely
%                  oriented dipole timecourse (6,T)
% T: number of timesamples
%--------------------------------------------------------------------------

dippos    = cfg.dippos; % dipole locations
headmodel = cfg.headmodel;
chantype  = cfg.channel;
unit      = cfg.unit;

% Extract dipole orientations in case they are not provided
%----------------------------------------------------------
if isfield(cfg,'ori')
    ori1 = cfg.ori(:,1);
    ori2 = cfg.ori(:,2);
else
    dipmom = cfg.dipmom;            % 6xT dipole moments used to extract orientation
    ori1   = mean(dipmom(1:3,:),2); % mean dipolmoment
    ori1   = ori1/norm(ori1);       % make unit vector
    ori2   = mean(dipmom(4:6,:),2); % mean dipolmoment
    ori2   = ori2/norm(ori2);       % make unit vector
end

% Extract channel type of interest
%---------------------------------
if contains(chantype,'meg')
    cfg         = [];
    cfg.channel = chantype;
    data        = ft_selectdata(cfg,data);
    sens        = data.grad;
elseif contains (chantype,'eeg')
    cfg         = [];
    cfg.channel = 'eeg*';
    data        = ft_selectdata(cfg,data);
    sens        = data.elec;
end

avg = data.avg;

% Extract Sphering information
%-----------------------------
if isfield(cfg,'noisecov')
    noisecov = cfg.noisecov;
else
    noisecov = [];
end

% Sphere data using the noise covariance matrix
%----------------------------------------------
if ~isempty(noisecov)
  disp('Sphering data with noise covariance.')
  [u, s]        = svd(noisecov);
  tol           = max(size(noisecov)) * eps(norm(s, inf));
  s             = diag(s);
  r1            = sum(s > tol) + 1;
  s(1:(r1 - 1)) = 1 ./ sqrt(s(1:(r1 - 1)));
  s(r1:end)     = 0;
  sphere        = diag(s) * u';
  % apply the sphering to the data
  avg           = sphere * avg;
  % apply the sphering as a pre-multiplication to the sensor definition
  montage          = [];
  montage.labelold = data.label;
  montage.labelnew = data.label;
  montage.tra      = sphere;
  sens             = ft_apply_montage(sens, montage, 'balancename', 'sphering');
end

% Computation of leadfield for predefined dipole positions
%---------------------------------------------------------
cfg                       = [];
cfg.grad                  = sens;  
cfg.channel               = chantype;                       
cfg.sourcemodel.pos       = dippos;
cfg.sourcemodel.unit      = unit;
cfg.headmodel             = headmodel;                
cfg.singleshell.batchsize = 5000;                      
lf                        = ft_prepare_leadfield(cfg);

% Apply orientation constraint to leadfield
%------------------------------------------
lf1 = lf.leadfield{1}*ori1;
lf2 = lf.leadfield{2}*ori2;
L   = [lf1,lf2];

% Solve overconstrained linear equation for scalar dipolmoments
%--------------------------------------------------------------
% avg = L*mom ,  avg (Nch,T), L (Nch,2), mom (2,T)
% mom = inv(L'*L)*L' * avg
% operator = inv(L'*L)*L';
operator      = (L'*L)\L';
scalar_dipmom = operator*avg; % scalar dipolmoments

end