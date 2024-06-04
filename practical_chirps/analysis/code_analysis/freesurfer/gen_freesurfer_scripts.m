close all
clear 
clc

%% Procedure
%--------------------------------------------------------------------------
% This script uses freesurfer and hence requires a Linux computer.
%
% 1.) Copy original BIDS data with anat folder in 
%     /home/till/Dokumente/MEG_T1/bids
%     e.g. /home/till/Dokumente/MEG_T1/bids/sub-01/anat
%
% 2.) Generate autorecon1 script. Is executed over all subjects.
%     This script is generated in /freesurfer1/scripts and can be executed
%     there. All paths in this file are absolte paths.
%
% 3.) After termination, process data with fieldtrip. Intermediate results
%     are stored in the respective folder.
%     /fieldtrip: input files for /freesurfer2 - tina-inspired freesurfer
%     pipeline
%     /freesurfer3: input files for fieldtrip inspired freesurfer pipeline
% 
% 4.) Genrate autorecon-all script and execute it for each subject
%     separetely in shell: bash freesurfer_reconall.sh
%     The scripts are generated in /freesurfer2/scripts and
%     /freesurfer3/scripts and can be executed in those folders cause all
%     paths in the scrips are absolute paths
%
% 5.) Process data with hcp workbench

%% Define subjects
%--------------------------------------------------------------------------
subjects = 1:4;
N_sub    = length(subjects);

datapath = '/home/till/Dokumente/MEG_T1/'; % on Linux
bidspath = '/home/till/Dokumente/MEG_T1/bids';

%% 1.) Copy data in bids folder
%--------------------------------------------------------------------------

%% 2.) Generate freesurfer scripts - autorecon1
%--------------------------------------------------------------------------

mkdir(fullfile(datapath,'freesurfer1')) % freesurfer autorecon1
mkdir(fullfile(datapath,'fieldtrip')) % intermediate fieldtrip results
mkdir(fullfile(datapath,'freesurfer2')) % freesurfer autorecon -all (tina inspired)
mkdir(fullfile(datapath,'freesurfer3')) % freesurfer autorecon -all (fieldtrip inspired)

% expert option for high resolution
exp = fopen(fullfile(datapath,'freesurfer1','expert.opts'),'wt');
fprintf(exp,'mris_inflate -n 15\n');
fclose(exp);

% freesurfer autorecon-1 script
mkdir(fullfile(datapath,'freesurfer1','scripts'))

fid_T1 = fopen(fullfile(datapath,'freesurfer1','scripts','freesurfer_autorecon1.sh'),'w');
fprintf(fid_T1,'export FREESURFER_HOME=/usr/local/freesurfer\n');
fprintf(fid_T1,'source $FREESURFER_HOME/SetUpFreeSurfer.sh\n\n');

for subidx = subjects

    subject = ['sub-',num2str(subidx,'%02d')];

    % Check if mri's exist
    if exist(fullfile(bidspath,subject,'anat'),'dir')

        series = dir([fullfile(bidspath,subject,'anat') filesep '*.nii']);

        fprintf(fid_T1,['export STUDY_DIR=' fullfile(bidspath,subject,'anat'),'\n']);
        fprintf(fid_T1,['export SUBJECTS_DIR=' fullfile(datapath,'freesurfer1'),'\n']);
        fprintf(fid_T1,'cd $SUBJECTS_DIR\n'); 

        % 1 T1 exist
        if length(series) == 1
            
             fprintf(fid_T1,['recon-all -autorecon1 -subjid /' subject ' -hires -i $STUDY_DIR/' [subject,'_run-01_T1w.nii'] ' -noskullstrip -expert  $SUBJECTS_DIR/expert.opts\n\n']); 
        % 2 T1 exist
        elseif length(series) == 2
            
             fprintf(fid_T1,['recon-all -autorecon1 -subjid /' subject ' -hires -i $STUDY_DIR/' [subject,'_run-01_T1w.nii'] ' -i $STUDY_DIR/' [subject,'_run-02_T1w.nii'] ' -expert $SUBJECTS_DIR/expert.opts\n\n']); 
        
        end

    end

end % for each subject

fclose(fid_T1);

%% 3.) After execution of script - fieldtrip processing
%------------------------------------------------------

for subidx = subjects

    subject   = ['sub-',num2str(subidx,'%02d')];

    % Check if mri's exist
    if exist(fullfile(bidspath,subject,'anat'),'dir')
        
        mri          = ft_read_mri(fullfile(datapath,'freesurfer1',subject,'/mri/T1.mgz'));
        mri.coordsys = 'ras'; 

        cfg            = [];
        cfg.resolution = 0.75;
        cfg.dim        = [320 320 320];
        cfg.yrange     = [-119.625,119.625]+20; % original size [-119.625,119.625] +20mm shift bachkwards for head centering
        mri_rs         = ft_volumereslice(cfg, mri);

        % check reslicing
        %----------------
%         cfg = [];
%         ft_sourceplot(cfg,mri)
%         ft_sourceplot(cfg,mri_rs)

        % for tina-inspired freesurfer pipeline
        mkdir(fullfile(datapath,'fieldtrip',subject))
        cfg           = [];
        cfg.filename  = fullfile(datapath,'fieldtrip',subject,'T1.mgz');
        cfg.filetype  = 'mgz';
        cfg.parameter = 'anatomy';
        ft_volumewrite(cfg, mri_rs);

        % for fieldtrip-inspired freesurfer pipeline
        cfg.filename  = fullfile(datapath,'freesurfer3',[subject,'.mgz']);
        cfg.filetype  = 'mgz';
        ft_volumewrite(cfg, mri_rs);

        % just to have it, can be loaded into fieldtrip in windows
        cfg.filename  = fullfile(datapath,'fieldtrip',subject,[subject,'.nii']);
        cfg.filetype  = 'nifti';
        ft_volumewrite(cfg, mri_rs);
        
        clear mri_rs mri
  
    end

end % loop over subjects

%% 4.) Generate freesurfer scripts 
%--------------------------------------------------------------------------

% One script per subject

% fieldtrip-inspired freesurfer pipeline
%---------------------------------------

mkdir(fullfile(datapath,'freesurfer3','scripts'))
for subidx = subjects

    subject   = ['sub-',num2str(subidx,'%02d')];

    % Check if mri's exist
    if exist(fullfile(bidspath,subject,'anat'),'dir')

        % Generate one script per participant
        fid = fopen(fullfile(datapath,'freesurfer3','scripts',['ft_freesurferscript_' matlab.lang.makeValidName(subject) '.sh']),'w');
        fprintf(fid,'export FREESURFER_HOME=/usr/local/freesurfer\n');
        fprintf(fid,'source $FREESURFER_HOME/SetUpFreeSurfer.sh\n');
        fprintf(fid,['export SUBJECTS_DIR=' datapath '/freesurfer3\n\n']);
    
        fprintf(fid,['mksubjdirs $SUBJECTS_DIR/' subject '\n']);
        fprintf(fid,'cd $SUBJECTS_DIR\n'); 
        fprintf(fid,['cp -f ' subject '.mgz $SUBJECTS_DIR/' subject '/mri/\n']); 
        fprintf(fid,['cd $SUBJECTS_DIR/' subject '/mri\n']); 
        fprintf(fid,['mri_convert -c -oc 0 0 0 ' subject '.mgz orig.mgz\n']);
        fprintf(fid,'cp orig.mgz orig/001.mgz\n\n');
    
        fprintf(fid,['recon-all -autorecon1 -subjid ' subject '\n']);
        fprintf(fid,['recon-all -autorecon2 -subjid ' subject '\n']);
        fprintf(fid,['recon-all -autorecon3 -subjid ' subject '\n\n']);
        fclose(fid);

    end

end

% bash freesurfer_reconall_fieldtrip.sh in subjectdir

% tina-inspired freesurfer pipeline
%----------------------------------

% expert option for high resolution
exp = fopen(fullfile(datapath,'freesurfer2','expert.opts'),'wt');
fprintf(exp,'mris_inflate -n 15\n');
fclose(exp);

mkdir(fullfile(datapath,'freesurfer2','scripts'))
for subidx = subjects

    subject   = ['sub-',num2str(subidx,'%02d')];

    % Check if mri's exist
    if exist(fullfile(bidspath,subject,'anat'),'dir')

        % Generate one script per participant
        fid = fopen(fullfile(datapath,'freesurfer2','scripts',['freesurfer_reconall_' matlab.lang.makeValidName(subject) '.sh']),'w');
        fprintf(fid,'export FREESURFER_HOME=/usr/local/freesurfer\n');
        fprintf(fid,'source $FREESURFER_HOME/SetUpFreeSurfer.sh\n');
        fprintf(fid,['export STUDY_DIR=' fullfile(datapath,'fieldtrip') '\n']);
        fprintf(fid,['export SUBJECTS_DIR=' fullfile(datapath,'freesurfer2') '\n\n']);

        fprintf(fid,'cd $SUBJECTS_DIR\n');   
%         fprintf(fid,['recon-all -all -subjid /' subject ' -hires -i $STUDY_DIR/' subject '/T1.mgz -expert $SUBJECTS_DIR/expert.opts\n\n']);
        fprintf(fid,['recon-all -all -subjid /' subject ' -hires -i $STUDY_DIR/' subject '/T1.mgz -expert $SUBJECTS_DIR/expert.opts -parallel -openmp 4\n\n']);
        fclose(fid);

    end

end

