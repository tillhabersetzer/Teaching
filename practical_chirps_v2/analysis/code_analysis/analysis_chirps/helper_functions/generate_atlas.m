function atlas = generate_atlas(atlasname,meshname,sourcemodel,subjectdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generates a surface based atlas from the freesurfer data
% it combines the atlases from both hemispheres
% subjectdata
% atlas name for left hemisphere
% - '.L.aparc.a2009s.4k_fs_LR.label.gii'
% - '.L.aparc.4k_fs_LR.label.gii'
% mesh name for left hemisphere
% - '.L.midthickness.4k_fs_LR.surf.gii'
% for 4k atlas or change to 8k
% sourcemodel is needed to get the order of the hemispheres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% paths
%------
path_label_l = fullfile(subjectdata.sourcemodel,'workbench',[subjectdata.subjectname,atlasname]);
path_mesh_l  = fullfile(subjectdata.sourcemodel,'workbench',[subjectdata.subjectname,meshname]);

% order
%------
label = sourcemodel.brainstructurelabel(sourcemodel.brainstructure(1));
if contains(label,'left','Ignorecase',true)
    atlas1 = ft_read_atlas({path_label_l,path_mesh_l});
    atlas2 = ft_read_atlas({strrep(path_label_l,'.L.','.R.'),strrep(path_mesh_l,'.L.','.R.')});
elseif contains(label,'right','Ignorecase',true)
    atlas1 = ft_read_atlas({strrep(path_label_l,'.L.','.R.'),strrep(path_mesh_l,'.L.','.R.')});
    atlas2 = ft_read_atlas({path_label_l,path_mesh_l});
else 
    error('Unsupported hemisphere order!')
end

% units needs to be changed in script where atlas is used
% atlas1 = ft_convert_units(atlas1,data_preprocessed_w.grad.unit);
% atlas2 = ft_convert_units(atlas2,data_preprocessed_w.grad.unit);

% combine both atlases from both hemispheres
%-------------------------------------------
L = length(atlas1.parcellationlabel);
P = length(atlas1.pos);

if contains(atlasname,{'.L.aparc.a2009s.4k_fs_LR.label.gii','.L.aparc.a2009s.8k_fs_LR.label.gii'})
    
    % check if atlas has known format
    if (~strcmp(atlas1.parcellationlabel(1),'???') || ~strcmp(atlas2.parcellationlabel(1),'???'))
        error('Cannot work with this atlas')
    end

    % change elements of the second hemissphere to add
    atlas.parcellationlabel = [atlas1.parcellationlabel;atlas2.parcellationlabel(2:end)];
    idx                     = ~ismember(atlas2.parcellation,1);
    parcel2change           = atlas2.parcellation;   
    % change parcellation numbers of atlas2
    % -1 because label '???' is removed on second atlas 
    parcel2change(idx)      = parcel2change(idx)+L-1; 
    atlas.parcellation      = [atlas1.parcellation;parcel2change];
    atlas.rgba              = [atlas1.rgba;atlas2.rgba(2:end,:)];
    
elseif contains(atlasname,{'.L.aparc.4k_fs_LR.label.gii','.L.aparc.8k_fs_LR.label.gii'})
    atlas.parcellationlabel = [atlas1.parcellationlabel;atlas2.parcellationlabel];
    atlas.parcellation      = [atlas1.parcellation;atlas2.parcellation+L];
    atlas.rgba              = [atlas1.rgba;atlas2.rgba];
end

atlas.pos               = [atlas1.pos;atlas2.pos];
atlas.tri               = [atlas1.tri;atlas2.tri+P];
atlas.unit              = atlas1.unit;

end