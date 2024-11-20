# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:36:33 2024

@author: Till Habersetzer
         Carl von Ossietzky University Oldenburg
         till.habersetzer@uol.de
"""

import os.path as op
import mne
from scipy.spatial.transform import Rotation as R
import numpy as np

# Compare transformations
#------------------------
subject = "sub-00"
task = 'transient_all'

# path to project (needs to be adjusted)
project_dir = r"C:\Users\tillhabersetzer\Nextcloud\Synchronisation\Projekte\GitHub\Oksima\rawdata\pilot\study"

#Load crosstalk compensation and fine calibration files
crosstalk_file = op.join(r'C:\Users\tillhabersetzer\Nextcloud\Synchronisation\Projekte\GitHub\Oksima\derivatives\SSS', 'ct_sparse.fif')
fine_cal_file = op.join(r'C:\Users\tillhabersetzer\Nextcloud\Synchronisation\Projekte\GitHub\Oksima\derivatives\SSS', 'sss_cal.dat')

#%% Compare device-to-head transformations for both sessions
#------------------------------------------------------------------------------
filename1 = subject + '_ses-01_task-' + task +'.fif'
data_path1 = op.join(project_dir,subject,'ses-01','meg',filename1)
raw1 = mne.io.read_raw_fif(data_path1, allow_maxshield=False, verbose=True)
trafo1 = raw1.info['dev_head_t']['trans']

filename2 = subject + '_ses-02_task-' + task +'.fif'
data_path2 = op.join(project_dir,subject,'ses-02','meg',filename2)
raw2 = mne.io.read_raw_fif(data_path2, allow_maxshield=False, verbose=True)
trafo2 = raw2.info['dev_head_t']['trans']


import sys
# Redirect stdout and stderr to a log file
log_file = open("console_output.log", "w", encoding="utf-8")
sys.stdout = log_file
sys.stderr = log_file


mne.set_log_file(None)     # Disable MNE's default file logging



# Use headposition of first recorded session as reference
destination   = data_path1
raw2_tsss_hpt = mne.preprocessing.maxwell_filter(raw2, cross_talk=crosstalk_file, 
                                                 calibration=fine_cal_file, st_duration=10, 
                                                 destination=destination, verbose=True)   
                                            