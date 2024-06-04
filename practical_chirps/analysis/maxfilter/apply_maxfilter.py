# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 13:06:49 2022

@author: tillhabersetzer
"""

#%% Settings
import os
import os.path as op
import mne
from mne.preprocessing import find_bad_channels_maxwell

subjects  = ['sub-01','sub-02']
tasks = ['clicks',
         'downchirps',
         'upchirps',
         'mmn',
         'emptyroom']

rootpath = os.path.join('M:',os.sep,'Blockpraktikum2022','analysis')

# Load crosstalk compensation and fine calibration files
crosstalk_file = os.path.join(rootpath,'derivatives','SSS', 'ct_sparse.fif')
fine_cal_file = os.path.join(rootpath,'derivatives','SSS', 'sss_cal.dat')

# Apply Oversampled Temporal Projection to reduce sensor noise before MaxFilter
OTP = 0
         
# Apply Headposition Transformation 
HPT = 1

#%% Apply Movement correction
MC = 0

#%% Headposition computations

if HPT: 
    for subject in subjects:
        
        figs_list= []
        captions_list = []
        dir2save =  os.path.join(rootpath,'derivatives',subject,'maxfilter')

        for task in tasks[0:4]:
            
            raw_fname = os.path.join(rootpath,'rawdata',subject,'meg',subject+'_task-' + task + '.fif')
            headpos_fname = os.path.join(rootpath,'derivatives',subject,'maxfilter',subject+'_task-' + task + '_head_pos.pos')
            
            # Head positions are computed if raw meg file exists and positions
            # havent been computed yet
            if op.isfile(raw_fname) and not op.isfile(headpos_fname):
                
                raw = mne.io.read_raw_fif(raw_fname, allow_maxshield=False, verbose=True)
            
                #%% Compute head position
                chpi_amplitudes = mne.chpi.compute_chpi_amplitudes(raw)
                chpi_locs = mne.chpi.compute_chpi_locs(raw.info, chpi_amplitudes)
                head_pos = mne.chpi.compute_head_pos(raw.info, chpi_locs, verbose=True)
                
                if head_pos.shape[0]>0: # cHPI active
                
                    #%% Save head position     
                    if not op.exists(dir2save):
                        os.makedirs(dir2save)
                        print("Directory '{}' created".format(dir2save))
        
                    mne.chpi.write_head_pos(headpos_fname, head_pos)
                    
                    captions_list.append(task)
                    figs_list.append(mne.viz.plot_head_positions(head_pos, mode='traces',show=False))
        
        # add to report if list is not empty
        if figs_list:
            #%% Add plots of the data to the HTML report
            report_fname = op.join(dir2save,subject+'-report.hdf5')
            report_html_fname = op.join(dir2save,subject+'-report.html')
            
            with mne.open_report(report_fname) as report:
                report.add_figure(
                figs_list,
                title='Extracting and visualizing subject head movement',
                caption=captions_list,
                replace=True
                )
            report.save(report_html_fname, overwrite=True,open_browser=False)  

#%% maxfilter processing

for subject in subjects:
    
    figs_list_before = []
    figs_list_after = []
    captions_list = []
    dir2save = os.path.join(rootpath,'derivatives',subject,'maxfilter')
    
    for task in tasks:
        
        #%% Load data
        raw_fname = os.path.join(rootpath,'rawdata',subject,'meg',subject+'_task-'+task+'.fif')
    
        if op.isfile(raw_fname):
            
            raw = mne.io.read_raw_fif(raw_fname, allow_maxshield=False, verbose=True)
        
            #%% Oversampled temporal projection
            if OTP:
                raw = mne.preprocessing.oversampled_temporal_projection(raw)
            
            #%% emptyroom 
            if task=='emptyroom':      
                destination = None
                head_pos = None
                st_duration = None
                coord_frame = "meg"
                    
            #%% recordings with subjects inside meg
            else:
                st_duration = 10
                coord_frame = 'head'
                
                #%% Head Position Transformation
                if HPT: 
                    # Use headposition of first recording as reference
                    ref_task = 'clicks'
                    destination = os.path.join(rootpath,'rawdata',subject,'meg',subject+'_task-'+task+'.fif')
                else:
                    destination = None
                    
                #%% Movement Correction  
                if MC: 
                    headpos_fname = os.path.join(rootpath,'derivatives',subject,'maxfilter',subject+'_task-' + task + '_head_pos.pos')
                    if op.isfile(headpos_fname):
                        head_pos = headpos_fname    
                else:
                    head_pos = None       
                
            #%% Detect bad channels
            raw.info['bads'] = []
            raw_check = raw.copy()
            auto_noisy_chs, auto_flat_chs = find_bad_channels_maxwell(
                raw_check, cross_talk=crosstalk_file, calibration=fine_cal_file,
                coord_frame=coord_frame, return_scores=False, verbose=True)
            print(auto_noisy_chs)  
            print(auto_flat_chs)  
            
            # Update list of bad channels
            bads = raw.info['bads'] + auto_noisy_chs + auto_flat_chs
            raw.info['bads'] = bads
            
            #%% Apply MaxFilter
            raw_tsss = mne.preprocessing.maxwell_filter(
                raw, cross_talk=crosstalk_file, calibration=fine_cal_file, 
                st_duration=st_duration, head_pos=head_pos, destination=destination,coord_frame=coord_frame, verbose=True)
                
            #%% Save data
            if not op.exists(dir2save):
                os.makedirs(dir2save)
                print("Directory '{}' created".format(dir2save))
                
            raw_tsss.save(os.path.join(dir2save,subject+'_task-'+task+'-raw_tsss.fif'),overwrite=True)

            #%% Add a plot of the data to the HTML report
            # report_fname = op.join(dir2save,subject+'-report.hdf5')
            # report_html_fname = op.join(dir2save,subject+'-report.html')
            # with mne.open_report(report_fname) as report:
            #     report.add_raw(raw=raw, title=f'Raw data before maxwell filter : {task}', psd=True, replace=True)
            #     report.add_raw(raw=raw_tsss, title=f'Raw data after maxwell filter: {task}', psd=True, replace=True)
            #     report.save(report_html_fname, overwrite=True,open_browser=False)
                
            figs_list_before.append(raw.plot_psd(show=False))
            figs_list_after.append(raw_tsss.plot_psd(show=False))
            captions_list.append(task)
        
    #%% Append plots to report
    report_fname = op.join(dir2save,subject+'-report.hdf5')
    report_html_fname = op.join(dir2save,subject+'-report.html')
    with mne.open_report(report_fname) as report:
        report.add_figure(
        figs_list_before,
        title='PSD before maxwell filtering',
        caption=captions_list,
        replace=True
        )
        report.add_figure(
        figs_list_after,
        title='PSD after maxwell filtering',
        caption=captions_list,
        replace=True
        )
    report.save(report_html_fname, overwrite=True,open_browser=False)  