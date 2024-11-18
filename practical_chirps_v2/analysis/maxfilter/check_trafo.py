# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:51:09 2024

@author: Till Habersetzer
         Carl von Ossietzky University Oldenburg
         till.habersetzer@uol.de
"""

import os.path as op
import mne
from scipy.spatial.transform import Rotation as R
import numpy as np

# Check transformations
#----------------------

subject = "sub-01"

tasks = [
    "ClickSensimetrics",
    "UpSensimetrics",
    "ClickTip300",
    "UpTip300",
]

# path to project (needs to be adjusted)
project_dir = r"C:\Users\tillhabersetzer\Nextcloud\Synchronisation\Projekte\GitHub\Teaching\practical_chirps_v2"


# Load data
dir2load1 = op.join(project_dir, "rawdata", subject, "meg")
dir2load2 = op.join(project_dir, "derivatives", subject, "maxfilter")

trafos1 = dict() # raw
trafos2 = dict() # maxfiltered
trafos_head1_to_head2 = dict() # head - head

for task in tasks:
    
    # Load both trafos and estimate translation and rotation
    #-------------------------------------------------------
      
    filename1 = subject + "_task-" + task + ".fif"
    filename2 = subject + "_task-" + task + "_proc-tsss" + "_meg.fif"
    
    data_path1 = op.join(dir2load1, filename1)
    raw1 = mne.io.read_raw_fif(data_path1, allow_maxshield=False, verbose=True)
    trafo1 = raw1.info['dev_head_t']['trans']
    trafos1[filename1] = trafo1
    
    data_path2 = op.join(dir2load2, filename2)
    raw2 = mne.io.read_raw_fif(data_path2, allow_maxshield=False, verbose=True)
    trafo2 = raw2.info['dev_head_t']['trans']
    trafos2[filename2] = trafo2
    
    # Compute the transformation from head1 to head2
    #-----------------------------------------------
    # 4x4 transformation matrices can be easily combined using matrix multiplication.
    # Since these matrices represent both rotation and translation in homogeneous coordinates, 
    # multiplying them correctly combines their transformations in sequence.
    T_head1_to_head2 = trafo2 @ np.linalg.inv(trafo1)  
    
    # Distance between maxfiltered and rawdata
    distance = np.linalg.norm(T_head1_to_head2[:3, 3])
    
    # Convert the rotation matrix to Euler angles (in radians)
    # Specify the order of rotations: "xyz" for rotations around x (pitch), y (yaw), z (roll)
    rotation = R.from_matrix(T_head1_to_head2[:3, :3])
    euler_angles = rotation.as_euler('xyz', degrees=True)
    # pitch, yaw, roll = euler_angles
    
    
    trafos_head1_to_head2[filename1] = {}
    trafos_head1_to_head2[filename1]['trafo'] = T_head1_to_head2
    trafos_head1_to_head2[filename1]['distance'] = distance
    trafos_head1_to_head2[filename1]['rotation'] = euler_angles
    
    
print('\nRawdata')
print('-------')
for filename, trafo in trafos1.items():  
    print(f"{filename}")
    print("-----------------------------------------------\n")
    print(f"{trafo}\n")
    
print('\nMaxfiltered data')  
print('----------------')
for filename, trafo in trafos2.items():  
    print(f"{filename}")
    print("-----------------------------------------------\n")
    print(f"{trafo}\n")

print('\nOverall transformation: head to head')  
print('------------------------------------')
for filename, trafo in trafos_head1_to_head2.items():  
    print(f"{filename}")
    print("-----------------------------------------------\n")
    print(f"{trafo['trafo']}\n")
    print(f"Translation distance: {trafo['distance']}\n")

    pitch, yaw, roll = trafo['rotation']
    print(f"Pitch (x-axis rotation): {pitch:.2f} degrees\n")
    print(f"Yaw (y-axis rotation): {yaw:.2f} degrees\n")
    print(f"Roll (z-axis rotation): {roll:.2f} degrees\n")
    
    
    

   
    
   
    
    









