% sub_02.m

% ensure that we don't mix up subjects
clear subjectdata

subjectdata.subjectnr  = '02';
subjectdata.age        = 'n/a';
subjectdata.sex        = 'm';
subjectdata.handedness = 'right';

subjectdata.comments   = ['Slight problems with microphone communication, ' ...
                          'MEG microphone cable was perhaps not plugged in correctly. ' ...
                          'Only the MEG was recorded.'];

% subject specific instructions
subjectdata.instructions = [];

% subjectdata.mrifile_info{1} = {'20211117'};
% subjectdata.mrifile_info{2} = {'2021-11-17T10:56:35Z'};
% subjectdata.mrifile_info{3} = {'DICOM_t1_mprage_sag_p2_iso75ausTRA_20211117105635_5.nii'};


