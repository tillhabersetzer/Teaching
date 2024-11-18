% sub_01.m

% ensure that we don't mix up subjects
clear subjectdata

subjectdata.subjectnr  = '01';
subjectdata.age        = 'n/a';
subjectdata.sex        = 'm';
subjectdata.handedness = 'right';

subjectdata.comments   = ['No cHPI during clicks, right earphone fits better '...
                          'than left (left: +20 dB gain necessary to reach hearing threshold). ' ...
                          'Only the MEG was recorded.'];

% subject specific instructions
subjectdata.instructions = [];

subjectdata.mrifile_info{1} = {};
subjectdata.mrifile_info{2} = {};
subjectdata.mrifile_info{3} = {};


