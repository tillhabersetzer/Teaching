# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:32:26 2024

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

print(f"{filename1}")
print("-----------------------------------------------\n")
print(f"{trafo1}\n")

print(f"{filename2}")
print("-----------------------------------------------\n")
print(f"{trafo2}\n")

# Confirm that both trafos are different
#---------------------------------------
# yes, same

# Estimate difference between trafos (rotation, translation)
#-----------------------------------------------------------
T_head1_to_head2 = trafo2 @ np.linalg.inv(trafo1)  

# Distance between two recordings
distance = np.linalg.norm(T_head1_to_head2[:3, 3])

# Convert the rotation matrix to Euler angles (in radians)
# Specify the order of rotations: "xyz" for rotations around x (pitch), y (yaw), z (roll)
# Euler angles define the angle of rotation around each respective axis
# extrinsic: rotations all refer to a fixed/global coordinate system xyz
# intrinsic: a rotation refers to the last rotated coordinate system 
# (starting with the first rotation that refers to the original/global coordinate system)
# use extrinsic rotation

rotation = R.from_matrix(T_head1_to_head2[:3, :3])
euler_angles = rotation.as_euler('xyz', degrees=True)
pitch, yaw, roll = euler_angles

print(f"Translation distance: {distance}\n")
print(f"Pitch (x-axis rotation): {pitch:.2f} degrees\n")
print(f"Yaw (y-axis rotation): {yaw:.2f} degrees\n")
print(f"Roll (z-axis rotation): {roll:.2f} degrees\n")

#%% Maxwell filtering (tSSS) with headposition transformation to first session
#------------------------------------------------------------------------------

# Use headposition of first recorded session as reference
destination   = data_path1
raw2_tsss_hpt = mne.preprocessing.maxwell_filter(raw2, cross_talk=crosstalk_file, 
                                                 calibration=fine_cal_file, st_duration=10, 
                                                 destination=destination, verbose=True)   
                                            
trafo2_tsss_hpt = raw2_tsss_hpt.info['dev_head_t']['trans']

print(f"{filename1}")
print("-----------------------------------------------\n")
print(f"{trafo1}\n")

print(f"{filename2}")
print("-----------------------------------------------\n")
print(f"{trafo2_tsss_hpt}\n")

# Confirm that trafos are the same now
#-------------------------------------
# yes, same

#%% Maxwell filtering (tSSS) with headposition transformation to default head position
#-------------------------------------------------------------------------------------

# compare with default head position
# destination=(0, 0, 0.04) would translate the bases as --trans default would 
# in MaxFilter™ (i.e., to the default head location).
raw1_tsss = mne.preprocessing.maxwell_filter(raw1, cross_talk=crosstalk_file, 
                                             calibration=fine_cal_file, st_duration=10,
                                             destination=(0, 0, 0.04),verbose=True)  
trafo1_tsss = raw1_tsss.info['dev_head_t']['trans']

raw2_tsss = mne.preprocessing.maxwell_filter(raw2, cross_talk=crosstalk_file, 
                                             calibration=fine_cal_file, st_duration=10,
                                             destination=(0, 0, 0.04),verbose=True)  
trafo2_tsss = raw2_tsss.info['dev_head_t']['trans']

print(f"{filename1}")
print("-----------------------------------------------\n")
print(f"{trafo1_tsss}\n")

print(f"{filename2}")
print("-----------------------------------------------\n")
print(f"{trafo2_tsss}\n")

# Estimate difference between trafos (rotation, translation)
#-----------------------------------------------------------
T_head1_to_head2 = trafo2_tsss @ np.linalg.inv(trafo1_tsss)  

# Distance between two recordings
distance = np.linalg.norm(T_head1_to_head2[:3, 3])

# Convert the rotation matrix to Euler angles (in radians)
# Specify the order of rotations: "xyz" for rotations around x (pitch), y (yaw), z (roll)
rotation = R.from_matrix(T_head1_to_head2[:3, :3])
euler_angles = rotation.as_euler('xyz', degrees=True)
pitch, yaw, roll = euler_angles

print(f"Translation distance: {distance}\n")
print(f"Pitch (x-axis rotation): {pitch:.2f} degrees\n")
print(f"Yaw (y-axis rotation): {yaw:.2f} degrees\n")
print(f"Roll (z-axis rotation): {roll:.2f} degrees\n")


# Confirm that both trafos are different
#---------------------------------------
# yes, same








