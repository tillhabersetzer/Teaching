% sub_04.m

% ensure that we don't mix up subjects
clear subjectdata

subjectdata.subjectnr  = '04';
subjectdata.age        = '25';
subjectdata.sex        = 'm';
subjectdata.handedness = 'right';

subjectdata.comments   = ['No empytroom recorded before measurement. ' ...
                          'MEG and 32 channel EEG were recorded.'];

% subject specific instructions
subjectdata.instructions = ["The participant was instructed to fixate a cross projected on a screen. "...
                            "He was instructed to keep his head and body still and blink as little as possible."];

% subjectdata.mrifile_info{1} = {'20200123','20200123'};
% subjectdata.mrifile_info{2} = {'2020-01-23T15:45:02Z','2020-01-23T15:51:18Z'};
% subjectdata.mrifile_info{3} = {'DICOM_t1_mprage_sag_p2_iso75ausTRA_20200123153215_5.nii','DICOM_t1_mprage_sag_p2_iso75ausTRA_20200123153215_6.nii'};

