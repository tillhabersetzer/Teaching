function mask = mne_generate_mask(roi,atlas_l,atlas_r,sourcemodel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a brain mask for a given sourcemodel based on atlases for each
% hemisphere.
% - hcp workbench atlases
% - mne sourcemodel
% - region of interests (roi) defined as a string or cell array of strings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do it for each hemisphere separately

% find correct parcellation number
label_atlas_l = find(contains(atlas_r.parcellationlabel,roi,'IgnoreCase',true));
label_atlas_r = find(contains(atlas_l.parcellationlabel,roi,'IgnoreCase',true));

% find correct surface points for given parcellation number
idx_atlas_l = ismember(atlas_l.parcellation,label_atlas_l);
idx_atlas_r = ismember(atlas_r.parcellation,label_atlas_r);

% check for left/right hemisphere order
label_hemisphere_l = find(contains(sourcemodel.brainstructurelabel,'Left','IgnoreCase',true));
label_hemisphere_r = find(contains(sourcemodel.brainstructurelabel,'Right','IgnoreCase',true));

% find correct hemisphere indices in sourcemodel
idx_hemisphere_l = find(ismember(sourcemodel.brainstructure,label_hemisphere_l));
idx_hemisphere_r = find(ismember(sourcemodel.brainstructure,label_hemisphere_r));

% create mask
mask = zeros(size(sourcemodel.inside));
mask(idx_hemisphere_l) = idx_atlas_l;
mask(idx_hemisphere_r) = idx_atlas_r;
mask = logical(mask);
end