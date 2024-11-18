# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 22:21:56 2022

@author: tillhabersetzer
"""

import os
import mne
from mne.preprocessing import find_bad_channels_maxwell

subjects  = ['sub-01','sub-02','sub-03','sub-04']

filenames = ['_task-noisepre_meg','_task-noisepost_meg']

# Apply Oversampled Temporal Projection to recude sensor noise before MaxFilter
OTP = 0 

path_root = os.path.join('M:',os.sep,'MEG_Lab_all','bids')

#### Load crosstalk compensation and fine calibration files
crosstalk_file = os.path.join(path_root,'derivatives','SSS', 'ct_sparse_160224_oldenburg.fif')
fine_cal_file = os.path.join(path_root,'derivatives','SSS', 'finecalibration_dc_160223_oldenburg.dat')

for subject in subjects:
    
    for filename in filenames:
        
        # check if file exists
        data_path = os.path.join(path_root,'sub-emptyroom','ses-sub'+subject[-2:],'meg','sub-emptyroom_ses-sub'+subject[-2:]+filename+'.fif')
        if os.path.isfile(data_path):

            # %% Load data
            
            raw = mne.io.read_raw_fif(data_path, allow_maxshield=False, verbose=True)
            
            # %% Oversampled temporal projection
            if OTP:
                raw = mne.preprocessing.oversampled_temporal_projection(raw)
            
            # %% Detect bad channels
    
            raw.info['bads'] = []
            raw_check = raw.copy()
            auto_noisy_chs, auto_flat_chs = find_bad_channels_maxwell(
                raw_check, cross_talk=crosstalk_file, calibration=fine_cal_file,
                return_scores=False, coord_frame="meg", verbose=True)
            print(auto_noisy_chs)  
            print(auto_flat_chs)  
    
            # Update list of bad channels
            bads = raw.info['bads'] + auto_noisy_chs + auto_flat_chs
            raw.info['bads'] = bads
            
            # %% Maxwell filtering (tSSS)
    
            raw_tsss = mne.preprocessing.maxwell_filter(
                raw, cross_talk=crosstalk_file, calibration=fine_cal_file, st_duration=10, coord_frame="meg", verbose=True)
            folder = 'tsss'
            
            # Save data
            directory = os.path.join(path_root,'derivatives',subject,folder)
            if not os.path.exists(directory):
                os.makedirs(directory)
                print("Directory '{}' created".format(folder))
    
            raw_tsss.save(os.path.join(directory,subject+filename+'-raw_tsss.fif'),overwrite=True)
    
            # Plot results
            # raw.pick(['meg']).plot(duration=10, highpass=1, lowpass=40, butterfly=True)
            # raw_tsss.pick(['meg']).plot(duration=10, highpass=1, lowpass=40, butterfly=True)
        