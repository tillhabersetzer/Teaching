function leadfield_sym = compute_sym_leadfield(cfg)
%--------------------------------------------------------------------------
% Till Habersetzer (10.02.2022)
% Within this script a symmetric leadfield model is computed. Therefore,
% those precomputed leadfields can be used in a symmetric dipolfit
% (symmetry in x-direction) to speed up the computation.
%
% Output:
%   leadfield_sym: symmetric leadfield model applicable to symmetric dipole
%                  fit.
%
% Input: 
%   cfg.sourcemodel: precomputed sourcemodel that will be made symmetric
%   cfg.headmodel:   precomputed headmodel
%   cfg.channel:     channels for leadfield computation
%   cfg.elec:        electrode definition (use either elec or grad)
%   cfg.grad:        gradiometer definition
%
% A few remarks:
% If you use the function before a call of ft_dipolfitting when using all
% meg (mag+grad) a sphering of the data is recommended. This sphering can 
% be done in ft_dipolfitting itself. 
% Unfortunately, using those precomputed leadfields with the sphering
% option doesnt work out, because in that case only the data is sphered but
% not the leadfields because those were precomputed here.
% This sphering option only makes sense if the leadfields are computed
% within the ft_dipolfitting function and sphered accordingly.
% That means, if you want to precompute the leadfields and use all meg
% channels together with sphering, you need to manually sphere the data
% before the call of this functions so that the leadfields are computed
% based on the updated magnetometer function. You can do this with:
%
% avg_sphered = avg_pooled;
% if sphering
%   noisecov = avg_noise.cov;
%   [u, sv] = svd(noisecov);
%   tol = max(size(noisecov)) * eps(norm(sv, inf));
%   sv = diag(sv);
%   r1 = sum(sv > tol) + 1;
%   sv(1:(r1 - 1)) = 1 ./ sqrt(sv(1:(r1 - 1)));
%   sv(r1:end)     = 0;
%   sphere = diag(sv) * u';
%   % apply the sphering to the data
%   avg_sphered.avg = sphere * avg_sphered.avg;
%   % apply the sphering as a pre-multiplication to the sensor definition
%   montage = [];
%   montage.labelold = avg_sphered.label;
%   montage.labelnew = avg_sphered.label;
%   montage.tra = sphere;
%   avg_sphered.grad = ft_apply_montage(avg_sphered.grad, montage, 'balancename', 'sphering');
% end
% .... compute leadfields based on avg_sphered.grad
% 
% Due to the fact that leadfields for meg can be computed very fast, there
% is no need to compute the leadfields beforehand. Therefore, the sphering
% option in ft_dipolfitting can be used without thinking.
%--------------------------------------------------------------------------

sourcemodel = cfg.sourcemodel;
headmodel   = cfg.headmodel;
channel     = cfg.channel;

% Differentiate between sensortypes
%----------------------------------
use_grad = false;
use_elec = false;
if isfield(cfg,'grad')
    use_grad = true;
    grad     = cfg.grad;
end
if isfield(cfg,'elec')
    use_elec = true;
    elec     = cfg.elec;
end

% Computation of symmetric sourcemodel
%-------------------------------------
cfg             = [];
if use_grad
    cfg.grad = grad;
end 
if use_elec
    cfg.elec = elec;
end
cfg.headmodel   = headmodel;
cfg.symmetry    = 'x'; % symmetry in x-direction
cfg.sourcemodel = sourcemodel;
sourcemodel_sym = ft_prepare_sourcemodel(cfg);
% -> positions have dim N-dipoles x 6 and will be split into two N-dipoles
% x 3 sourcemodels for leadfield computation. The leadfields will be merged
% afterwards.

% Instantiate sourcemodels
%-------------------------
sourcemodel1     = sourcemodel;
sourcemodel1.pos = sourcemodel_sym.pos(:,1:3);
sourcemodel2     = sourcemodel;
sourcemodel2.pos = sourcemodel_sym.pos(:,4:6);

% Computation of leadfiels
%-------------------------
cfg                       = [];
if use_grad
    cfg.grad = grad;
end 
if use_elec
    cfg.elec = elec;
end
cfg.sourcemodel           = sourcemodel1;
cfg.headmodel             = headmodel;
cfg.channel               = channel;
cfg.singleshell.batchsize = 'all'; % number of dipoles for which the leadfield is computed in a single call
leadfield1                = ft_prepare_leadfield(cfg); 
cfg.sourcemodel           = sourcemodel2;
leadfield2                = ft_prepare_leadfield(cfg); 

% Merge leadfields
%-----------------
leadfield_sym     = leadfield1;
leadfield_sym.pos = [leadfield1.pos,leadfield2.pos];

for i = 1:length(leadfield_sym.inside)
    if leadfield_sym.inside(i)
        leadfield_sym.leadfield{i} = [leadfield1.leadfield{i},leadfield2.leadfield{i}];
    end
end

end