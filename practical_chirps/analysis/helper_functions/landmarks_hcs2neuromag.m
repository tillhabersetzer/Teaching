function [pos, label] = landmarks_hcs2neuromag(cfg)
%--------------------------------------------------------------------------
% Till Habersetzer (07.02.22)
%
% This script is designed to support the coregistration between MRI and MEG
% for BIDS. Therfore, anatomical landmarks (NAS, LPA, RPA) need to be
% determined in the head coordinate system (hcs) and in the MRI coordinate 
% system (voxel space). 
% The anatomical landmarks in hcs are measured during the digitization and 
% stored in the meg file and can be retrieved in the headshape object.
% The anatomical landmarks in the mri have to be determined by manual
% inspection. Typically, the first coregistration is refined afterwards by
% including the headshape. This further refinement leads to a missalignment
% between anatomical landmarks in the mri and hcs. Therefore, following the
% mne approach, the digitized landmarks are transformed into the mri voxel
% coordinate system after the coregistration and replace the previously 
% manually determined landmarks in the mri. With this approach, digitzed 
% landmarks and landmarks in the mri coincide.
% More detailed information can be found here: 
% https://groups.google.com/g/bids-discussion/c/BeyUeuNGl7I
%
% To be bids compliant, the coregistration needs to be stored in two
% places:
% - in _coordsystem.json: contains digitzed landmarks in hcs 
%                         same for each MEG session, stored in meg folder
% - in  _t1w.json: contains voxel coordinates of anatomical landmarks,
%                  stored in anat folder
%
% Input:
%   Shape and transform must be specified in same units.
%   cfg      : Configuration structure
%   cfg.shape: Headshape object calculated via ft_read_headshape, contains
%              labels and positions in hcs of anatomical landmarks.
%   cfg.transform: Homogeneous transformation matrix of coregistered mri.
%                  Transforms voxel coordinates of mri into neuromag
%                  coordinate system (hcs), (4x4 matrix)
% Output:
%   pos: Positions of anatomical landmarks in voxel coordinates
%   label: Labels of anatomical landmarks
%
% Application of coregistration:
% Use: ft_volumerealign
% 
% cfg              = [];
% cfg.coordsys     = 'neuromag';
% cfg.fiducial.nas = pos(2,:); % position of landmarks, order may change!
% cfg.fiducial.lpa = pos(1,:);
% cfg.fiducial.rpa = pos(3,:);
% mri_realigned    = ft_volumerealign(cfg,mri_orig);
% 
% for a fast check:
% cfg                       = [];
% cfg.method                = 'headshape';
% cfg.headshape.headshape   = shape;
% cfg.coordsys              = 'neuromag';
% cfg.headshape.interactive = 'yes';
% cfg.headshape.icp         = 'no'; 
% mri_realigned2            = ft_volumerealign(cfg, mri_realigned);
% 
%--------------------------------------------------------------------------

% Rename data
%------------
headshape         = cfg.shape;
fid               = headshape.fid; % extract anatomical landmarks, fiducials
traf_vox2neuromag = cfg.transform;
label             = fid.label;

% Apply homogeneous coordinae transformation
%--------------------------------------------------------------------------
% multiple ways possible

% 1.) Manual transformation
%--------------------------
% reorder coordinates, positions are stored row-wise originally and are 
% changed into homogeneous coordinates.
pos_neuromag = [fid.pos,[1;1;1]]';

% pos_vox = inv(traf_vox2neuromag)*pos_neuromag;
pos_vox = traf_vox2neuromag\pos_neuromag;
pos_vox = pos_vox(1:3,:)'; % reorder
pos     = pos_vox;
% pos     = round(pos_vox); % round to intergers for voxel coordinates

% 2.) Use fieldtrip function
%---------------------------
% apply homogenous coordinate transformation to shape
% shape_vox = ft_transform_geometry(inv(traf_vox2neuromag), shape);
% pos       = round(shape_vox.fid.pos);
% label     = shape_vox.fid.label;

end