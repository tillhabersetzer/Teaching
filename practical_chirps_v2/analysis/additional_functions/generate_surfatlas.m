function atlas = generate_surfatlas(atlasname,sourcemodel,rootpath,subject)
%--------------------------------------------------------------------------
% Till Habersetzer (05.02.2022)
% Generates a surface based atlas from the freesurfer and hcp workbench
% preprocessed data. It combines the atlases from both hemispheres.
% 
% atlas name for left hemisphere can be:
% - '.L.aparc.a2009s.4k_fs_LR.label.gii'
% - '.L.aparc.4k_fs_LR.label.gii'
% mesh name for left hemisphere:
% - '.L.midthickness.4k_fs_LR.surf.gii'
% for 4k atlas or change to 8k
% sourcemodel is needed to get the order of the hemispheres
%--------------------------------------------------------------------------

% Define paths to acces mesh and atlas
%-------------------------------------
if contains(atlasname,'4k')
    res = '4k';
elseif contains(atlasname,'8k')
    res = '8k';
else
    error('Atlas not supported.')
end
meshname = ['.L.midthickness.',res,'_fs_LR.surf.gii'];

% Check order of hemispheres so that order is identical to sourcemodel
%---------------------------------------------------------------------
label = sourcemodel.brainstructurelabel(sourcemodel.brainstructure(1));
path_label_l = fullfile(rootpath,'derivatives',subject,'freesurfer','workbench',[subject,atlasname]);
path_mesh_l  = fullfile(rootpath,'derivatives',subject,'freesurfer','workbench',[subject,meshname]);

% Load atlas for left and right hemisphere
%-----------------------------------------
if contains(label,'left','Ignorecase',true)
    atlas1 = ft_read_atlas({path_label_l,path_mesh_l});
    atlas2 = ft_read_atlas({strrep(path_label_l,'.L.','.R.'),strrep(path_mesh_l,'.L.','.R.')});
elseif contains(label,'right','Ignorecase',true)
    atlas1 = ft_read_atlas({strrep(path_label_l,'.L.','.R.'),strrep(path_mesh_l,'.L.','.R.')});
    atlas2 = ft_read_atlas({path_label_l,path_mesh_l});
else 
    error('Unsupported hemisphere order!')
end

% Combine atlases from both hemispheres
%--------------------------------------------------------------------------
L = length(atlas1.parcellationlabel); 
P = length(atlas1.pos);

% Differentiate between two availabe atlases
%-------------------------------------------
if contains(atlasname,{'.L.aparc.a2009s.4k_fs_LR.label.gii','.L.aparc.a2009s.8k_fs_LR.label.gii'})
    
    % Check if atlas has known format
    if (~strcmp(atlas1.parcellationlabel(1),'???') || ~strcmp(atlas2.parcellationlabel(1),'???'))
        error('Cannot work with this atlas.')
    end

    % Change parcellation of 2nd hemissphere
    atlas.parcellationlabel = [atlas1.parcellationlabel;atlas2.parcellationlabel(2:end)];
    idx                     = ~ismember(atlas2.parcellation,1); % change all elements except elements marked with 1
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