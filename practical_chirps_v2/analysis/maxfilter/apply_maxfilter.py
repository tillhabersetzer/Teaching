# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 18:34:16 2024

@author: Till Habersetzer
         Carl von Ossietzky University Oldenburg
         till.habersetzer@uol.de

Application of Maxwell filtering
https://mne.tools/stable/auto_tutorials/preprocessing/60_maxwell_filtering_sss.html

Optional:
---------
- Application of oversampled temporal projection (otp) (takes too long...)
  Denoising algorithm
  https://mne.tools/stable/auto_examples/preprocessing/otp.html
- Computation and correction of head movements
  https://mne.tools/dev/auto_tutorials/preprocessing/59_head_positions.html
- Transformation to common head positions between runs. That means same
  head-dev-trafo for all runs
"""

# %% Settings
import os
from pathlib import Path
import mne
from mne.preprocessing import find_bad_channels_maxwell

subjects = ["sub-01", "sub-02", "sub-03"]

fnames = [
    "ClickSensimetrics",
    "UpSensimetrics",
    "ClickTip300",
    "UpTip300",
    "emptyroom",
]

# path to project (needs to be adjusted)
project_dir = Path(
    r"C:\Users\tillhabersetzer\Nextcloud\Synchronisation\Projekte\GitHub\Teaching\practical_chirps_v2"
)  # / -> makes it an absolute path

# Load crosstalk compensation and fine calibration files
crosstalk_file = Path(project_dir, "derivatives", "SSS", "ct_sparse.fif")
fine_cal_file = Path(project_dir, "derivatives", "SSS", "sss_cal.dat")

# Apply Oversampled Temporal Projection to reduce sensor noise before MaxFilter
OTP = False

# Apply Headposition Transformation
HPT = True
ref_fname = "UpTip300"

# Apply and compute movement correction
MC = False

badchan2add = ["MEG1813"]

# %% Headposition computations
# ------------------------------------------------------------------------------
# Compute headposition transformations

if MC:
    for subject in subjects:
        figs_list = []
        captions_list = []

        dir2save = Path(project_dir, "derivatives", subject, "maxfilter")
        rawdata_dir = Path(project_dir, "rawdata", subject, "meg")

        # Save head position
        # -------------------
        if not dir2save.is_dir():
            os.makedirs(dir2save)
            print("Directory '{}' created".format(dir2save))

        # use reduced set for head movement computations
        fnames_reduced = [fname for fname in fnames if not "emptyroom" in fname]
        for fname in fnames_reduced:
            fname_prefix = subject + "_task-" + fname
            raw_fname = Path(rawdata_dir, fname_prefix + ".fif")
            headpos_fname = Path(dir2save, fname_prefix + "_raw.pos")

            # Head positions are computed if raw meg file exists and positions
            # havent been computed yet
            if raw_fname.is_file() and not headpos_fname.is_file():
                raw = mne.io.read_raw_fif(
                    raw_fname, allow_maxshield=False, verbose=True
                )

                # Compute head position
                # ----------------------
                chpi_amplitudes = mne.chpi.compute_chpi_amplitudes(raw)
                chpi_locs = mne.chpi.compute_chpi_locs(raw.info, chpi_amplitudes)
                head_pos = mne.chpi.compute_head_pos(raw.info, chpi_locs, verbose=True)

                # Check if continuous head position is available
                # -----------------------------------------------
                if head_pos.shape[0] > 0:  # cHPI active
                    mne.chpi.write_head_pos(headpos_fname, head_pos)

                    captions_list.append(fname)
                    figs_list.append(
                        mne.viz.plot_head_positions(head_pos, mode="traces", show=False)
                    )

        # Add to report if list is not empty
        # -----------------------------------
        if figs_list:
            # Add plots of the data to the HTML report
            report_fname = Path(dir2save, subject + "-report.hdf5")
            report_html_fname = Path(dir2save, subject + "-report.html")

            with mne.open_report(report_fname) as report:
                report.add_figure(
                    figs_list,
                    title="Extracting and visualizing subject head movement",
                    caption=captions_list,
                    replace=True,
                )
            report.save(report_html_fname, overwrite=True, open_browser=False)

# %% maxfilter processing
# ------------------------------------------------------------------------------

for subject in subjects:
    figs_list_before = []
    figs_list_after = []
    captions_list = []

    dir2save = Path(project_dir, "derivatives", subject, "maxfilter")
    rawdata_dir = Path(project_dir, "rawdata", subject, "meg")

    # Save head position
    # -------------------
    if not dir2save.is_dir():
        os.makedirs(dir2save)
        print("Directory '{}' created".format(dir2save))

    for fname in fnames:
        # Load data
        # ----------
        fname_prefix = subject + "_task-" + fname
        raw_fname = Path(rawdata_dir, fname_prefix + ".fif")

        if raw_fname.is_file():
            raw = mne.io.read_raw_fif(raw_fname, allow_maxshield=False, verbose=True)

            # Oversampled temporal projection
            # --------------------------------
            if OTP:
                raw = mne.preprocessing.oversampled_temporal_projection(raw)

            # emptyroom
            # ----------
            if "emptyroom" in fname:
                destination = None
                head_pos = None
                st_duration = None
                coord_frame = "meg"
                proc = ""

            # recordings with subjects inside meg
            # ------------------------------------
            else:
                st_duration = 10
                coord_frame = "head"

                # Head Position Transformation
                # -----------------------------
                if HPT:
                    # Use headposition of first recording as reference
                    destination = Path(
                        rawdata_dir,
                        subject + "_task-" + ref_fname + ".fif",
                    )
                else:
                    destination = None

                # Movement Correction
                # --------------------
                headpos_fname = Path(dir2save, fname_prefix + "_raw.pos")
                if MC and headpos_fname.is_file():
                    head_pos = headpos_fname
                    proc = "_proc-tsss+mc"
                else:
                    head_pos = None
                    proc = "_proc-tsss"

            # Detect bad channels
            # --------------------
            raw.info["bads"] = []
            raw_check = raw.copy()
            auto_noisy_chs, auto_flat_chs = find_bad_channels_maxwell(
                raw_check,
                cross_talk=crosstalk_file,
                calibration=fine_cal_file,
                coord_frame=coord_frame,
                return_scores=False,
                verbose=True,
            )
            print(auto_noisy_chs)
            print(auto_flat_chs)

            # Update list of bad channels
            bads = raw.info["bads"] + auto_noisy_chs + auto_flat_chs
            if badchan2add not in bads:
                bads += badchan2add
                
            raw.info["bads"] = bads
          
            # Apply MaxFilter
            # ----------------
            raw_tsss = mne.preprocessing.maxwell_filter(
                raw,
                cross_talk=crosstalk_file,
                calibration=fine_cal_file,
                st_duration=st_duration,
                head_pos=head_pos,
                destination=destination,
                coord_frame=coord_frame,
                verbose=True,
            )

            # Save data
            # ----------
            raw_tsss.save(
                Path(dir2save, fname_prefix + proc + "_meg.fif"), overwrite=True
            )

            # Add a plot of the data to the HTML report
            # ------------------------------------------
            # report_fname = op.join(dir2save,subject+'-report.hdf5')
            # report_html_fname = op.join(dir2save,subject+'-report.html')
            # with mne.open_report(report_fname) as report:
            #     report.add_raw(raw=raw, title=f'Raw data before maxwell filter : {task}', psd=True, replace=True)
            #     report.add_raw(raw=raw_tsss, title=f'Raw data after maxwell filter: {task}', psd=True, replace=True)
            #     report.save(report_html_fname, overwrite=True,open_browser=False)

            figs_list_before.append(raw.compute_psd().plot(show=False, xscale="log"))
            figs_list_after.append(
                raw_tsss.compute_psd().plot(show=False, xscale="log")
            )
            captions_list.append(fname)

    # Append plots to report
    # -----------------------
    report_fname = Path(dir2save, subject + "-report.hdf5")
    report_html_fname = Path(dir2save, subject + "-report.html")
    with mne.open_report(report_fname) as report:
        report.add_figure(
            figs_list_before,
            title="PSD before maxwell filtering",
            caption=captions_list,
            replace=True,
        )
        report.add_figure(
            figs_list_after,
            title="PSD after maxwell filtering",
            caption=captions_list,
            replace=True,
        )
    report.save(report_html_fname, overwrite=True, open_browser=False)
