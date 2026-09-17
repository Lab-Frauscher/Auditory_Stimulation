"""
thalamus_analysis.py

Description:
==================
This script processes thalamic stereo-EEG (SEEG) recordings in European Data Format (EDF) 
to quantify spectral power dynamics across auditory stimulation events. It evaluates 
bipolar differential signals constructed from specified thalamic contact pairs and extracts 
matched 3-second pre-stimulus (baseline) and post-stimulus/arousal epoch windows. 
Using zero-phase bandpass Butterworth filters, it computes delta-band (0.5–4.0 Hz) to 
higher-frequency (4.5–40.0 Hz) power ratios to characterize subcortical spectral shifts.

Key Processing Steps:
-------------------
1. Loading event annotations (arousals/controls) and thalamic channel mapping configuration.
2. Iterating through patient cohort EDF files and opening streams via PyEDFlib.
3. Resolving target bipolar electrode contacts and mapping channel index locations.
4. Extracting signal sampling frequency and computing UTC epoch origin timestamps (with explicit timezone handling).
5. Constructing conditional pre-stimulus (baseline) and post-stimulus (evoked/control) analysis windows.
6. Validating sample boundaries to prevent out-of-bounds read requests in C-backend.
7. Extracting raw channel epoch slices and calculating bipolar differential signals.
8. Bandpass filtering via zero-phase Butterworth filters for Delta (0.5-4.0 Hz) and 
   Theta/Alpha/Beta/Gamma (4.5-40.0 Hz) bands.
9. Feature Extraction: Computing median instantaneous power ratios (Delta / Higher-bands) 
   with zero-division guards for baseline and post-stimulus epochs.
10. Structuring output metrics into a consolidated pandas DataFrame for scientific reporting.

Dependencies:
-------------------
- pyedflib
- scipy
- pandas
- numpy
"""

from datetime import timezone
import numpy as np
import pandas as pd
import pyedflib
from scipy.signal import butter, filtfilt

# --- 1. DATA INITIALIZATION AND ANNOTATIONS LOADING ---
# Load arousal annotations and add incremental event numbering
arousals = pd.read_excel('arousal_annotations.xlsx')
arousals['Event number'] = range(1, len(arousals) + 1)

# Load target thalamic bipolar channel pairs mapping
channels = pd.read_excel('thalamus_channels.xlsx', sheet_name='channels_bipolar')

# Initialize target DataFrame for aggregated power ratio metrics
power_ratio = pd.DataFrame([])

# --- 2. THALAMIC CHANNEL PAIRS PROCESSING PIPELINE ---
for i, thal_chan_row in channels.iterrows():
    patient = thal_chan_row['patient']
    id_pacient = thal_chan_row['ID']
    thal_chan1 = str(thal_chan_row['channel1']).strip()
    thal_chan2 = str(thal_chan_row['channel2']).strip()

    # Construct target EDF file path
    file_name = f'demo_data/auditory_stimulation_{patient}_demo.edf'
    
    # Filter event annotations corresponding to current patient
    arousal_p = arousals.loc[arousals['ID'] == patient]
    
    # --- 3. RECORDING LOAD AND CHANNEL INDEX MAPPING ---
    try:
        f = pyedflib.EdfReader(file_name)
    except (OSError, IOError) as e:
        print(f"Error opening EDF file {file_name}: {e}. Skipping configuration.")
        continue

    try:
        unit_conversion = 1.0  # Unit scale factor (default: volts/microvolts)
        
        # Extract signal labels and match exact channel indices
        ch_labels = [label.strip() for label in f.getSignalLabels()]
        if thal_chan1 not in ch_labels or thal_chan2 not in ch_labels:
            print(f"Patient {patient}: Channel pair ({thal_chan1}, {thal_chan2}) not found. Skipping.")
            continue

        idx1 = ch_labels.index(thal_chan1)
        idx2 = ch_labels.index(thal_chan2)

        # Extract sampling frequency, total channel samples, and explicit UTC origin timestamp
        fvz = f.getSampleFrequency(idx1)
        total_samples = f.getNSamples()[idx1]
        
        start_datetime = f.getStartdatetime()
        if start_datetime is not None:
            # Force UTC interpretation on naive datetime objects to match annotation timestamps
            if start_datetime.tzinfo is None:
                start_datetime = start_datetime.replace(tzinfo=timezone.utc)
            time_start = start_datetime.timestamp()
        else:
            time_start = 0
        
        # --- 4. PRE/POST-STIMULUS EPOCH EXTRACTION & FEATURE COMPUTATION ---
        for idx, row in arousal_p.iterrows():
            a_type = row['type']
            t_stim = row['utc_stim_start']
            t_arou = row['utc_arousal_start']
            ev_num = row['Event number']
            stim_intensity = row.get('stim_intensity', np.nan)
            
            sleep_stage_pre = row.get('stage', np.nan)
            sleep_stage_post = row.get('stage_post', np.nan)
            
            # Define baseline (pre) and event (post) analysis windows based on event condition
            if a_type == 'evoked':
                pre_start, pre_end = t_stim - 3, t_stim
                post_start, post_end = t_arou, t_arou + 3
            elif a_type == 'control':
                pre_start, pre_end = t_stim - 3, t_stim
                post_start, post_end = t_stim + 1.076, t_stim + 1.076 + 3
            else: 
                continue
            
            # Convert absolute UTC timestamp to relative recording time (seconds)
            post_tmin, post_tmax = post_start - time_start, post_end - time_start
            
            # Convert relative time window (s) to PyEDFlib sample indices
            post_start_sample = int(np.round(post_tmin * fvz))
            post_n_samples = int(np.round((post_tmax - post_tmin) * fvz))
            
            # Validate sample index bounds against recording length
            if post_start_sample < 0 or (post_start_sample + post_n_samples) > total_samples:
                print(f"Patient {patient}, Event {ev_num}: Post-stimulus window out of file bounds. Skipping.")
                continue

            # Extract raw signal slices for post-stimulus / arousal epoch
            data1 = f.readSignal(idx1, post_start_sample, post_n_samples)
            data1 = -1 * data1 * unit_conversion
            
            data2 = f.readSignal(idx2, post_start_sample, post_n_samples)
            data2 = -1 * data2 * unit_conversion

            # Construct bipolar differential signal
            signal = (data1 - data2).ravel()

            # Filter 1: Delta frequency band (0.5 - 4.0 Hz, 2nd order Butterworth, zero-phase)
            b1, a1 = butter(2, [0.5 / (fvz / 2), 4.0 / (fvz / 2)], 'bandpass')
            signal_d = filtfilt(b1, a1, signal)
            
            # Filter 2: Theta/Alpha/Beta/Gamma frequency band (4.5 - 40.0 Hz, 3rd order Butterworth, zero-phase)
            b2, a2 = butter(3, [4.5 / (fvz / 2), 40.0 / (fvz / 2)], 'bandpass')
            signal_tab = filtfilt(b2, a2, signal)
                    
            # Feature Computation: Ratio of median squared power (Delta / Higher Bands) with zero protection
            denom_post = np.median(signal_tab ** 2)
            result_powerratio = (np.median(signal_d ** 2) / denom_post) if denom_post > 0 else np.nan
            
            # --- BASELINE (PRE-STIMULUS) EPOCH EXTRACTION ---
            pre_tmin, pre_tmax = pre_start - time_start, pre_end - time_start
            
            pre_start_sample = int(np.round(pre_tmin * fvz))
            pre_n_samples = int(np.round((pre_tmax - pre_tmin) * fvz))
            
            # Validate sample index bounds for baseline epoch
            if pre_start_sample < 0 or (pre_start_sample + pre_n_samples) > total_samples:
                print(f"Patient {patient}, Event {ev_num}: Pre-stimulus window out of file bounds. Skipping.")
                continue

            # Extract raw signal slices for pre-stimulus baseline epoch
            data_base1 = f.readSignal(idx1, pre_start_sample, pre_n_samples)
            data_base1 = -1 * data_base1 * unit_conversion
            
            data_base2 = f.readSignal(idx2, pre_start_sample, pre_n_samples)
            data_base2 = -1 * data_base2 * unit_conversion

            # Construct bipolar differential signal for baseline
            signal_base = (data_base1 - data_base2).ravel()

            # Filter baseline signals across defined frequency bands
            signal_d_base = filtfilt(b1, a1, signal_base)
            signal_tab_base = filtfilt(b2, a2, signal_base)
                    
            # Feature Computation: Baseline spectral power ratio with zero protection
            denom_base = np.median(signal_tab_base ** 2)
            result_powerratio_base = (np.median(signal_d_base ** 2) / denom_base) if denom_base > 0 else np.nan

            # Aggregate trial metric row
            power_ratio = pd.concat([power_ratio, pd.DataFrame({
                'Patient ID': patient,
                'Event type': a_type,
                'Event number': ev_num,
                'Stimulus name': stim_intensity,
                'Pre-stimulus Sleep Stage': sleep_stage_pre,
                'Post-stimulus Sleep Stage': sleep_stage_post,
                'SEEG channel': [f"{thal_chan1}-{thal_chan2}"],
                'Pre-stimulus spectral power ratio': result_powerratio_base,
                'Arousal (/event) spectral power ratio': result_powerratio
            })], ignore_index=True)

    finally:
        f.close()

# --- 5. EXPORT / STORAGE (OPTIONAL) ---
# power_ratio.to_excel('thalamus_power_ratio_results.xlsx', index=False)
