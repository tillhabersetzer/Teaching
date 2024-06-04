% sub_03.m

% ensure that we don't mix up subjects
clear subjectdata

subjectdata.subjectnr  = '03';
subjectdata.age        = 'n/a';
subjectdata.sex        = 'm';
subjectdata.handedness = 'right';

subjectdata.comments   = ['Probably more noise on right hemisphere sensors, ' ...
                          'The participant cleared his throat occasionally. '...
                          'Only the MEG was recorded.'];

% subject specific instructions
subjectdata.instructions = [];

% subjectdata.mrifile_info{1} = {'20211117'};
% subjectdata.mrifile_info{2} = {'2021-11-17T10:45:37Z'};
% subjectdata.mrifile_info{3} = {'DICOM_t1_mprage_sag_p2_iso75ausTRA_20211117104537_5.nii'};