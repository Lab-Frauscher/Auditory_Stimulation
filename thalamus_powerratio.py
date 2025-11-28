"""
thalamus_powerratio.py

Description:
==================

This script implements a signal processing and feature extraction pipeline designed
to investigate thalamocortical activation. The analysis is performed on a subset
of patients with SEEG electrodes implanted within the thalamus. The core feature
calculated is the spectral power ratio of the thalamic intracranial EEG (iEEG) 
signal, specifically δ/(θ+α+β), which covers the range of [0.5–4 Hz] / [4.5–40 Hz].
This ratio is computed over the 3-second arousal epoch and then normalized to 
a baseline period defined before the stimulus onset.

Key Processing Steps:
-------------------
1. Loading arousal times and thalamic channel names 
2. Accessing iEEG signal data (.mat files)
3. Extracting the signal for the 3-second arousal epoch and the 3-second baseline epoch
4. Signal filtering to the Delta frequency band and the Theta/Alpha/Beta band.
5. Power ratio calculation for both the arousal and baseline periods.
6. Relative change of power ratio divided by the baseline power ratio to determine
the relative change
7. Aggregating results to find the average change per arousal across all thalamic channels

Dependencies:
-------------------
- numpy
- pandas
- scipy.signal
- h5py
"""

import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt
import h5py

# Loading arousal annotations and thalamic channels
arousals = pd.read_excel('thalamus_channels.xlsx',sheet_name='arousal')
channels = pd.read_excel('thalamus_channels.xlsx',sheet_name='channels_bipolar')

# Set window length for 3 seconds
arousal_duration_seconds = 3

power_ratio = pd.DataFrame([])

# Iterate over thalamus channels
for i,thal_chan_row in channels.iterrows():
    patient = thal_chan_row['patient']
    thal_chan = thal_chan_row['channel']
    arousal_p = arousals.loc[arousals['patient']==patient]

    # File name
    file_path = 'demo_data/'+patient+'.mat'

    # The list of variables to extract
    variables_to_load = ['bp_ch_names', 'bp_data', 'fs', 'n_bp']
    
    # Process .mat file for Python use
    with h5py.File(file_path, 'r') as f:       
        for var_name in variables_to_load:
            data = f[var_name][:] 
            
            # For string data
            if var_name == 'bp_ch_names':
                ch_names = []
                for ref in data.flatten():
                    if ref:
                        name_bytes = f[ref][()].tobytes()
                        cleaned_bytes = name_bytes.replace(b'\x00', b'')
                        ch_names.append(cleaned_bytes.decode('utf-8'))
            # For numerical data
            else:
                if var_name == 'bp_data':
                    bp_data = np.squeeze(data)
                
                elif var_name == 'fs':
                    fs = np.squeeze(data)
    
    # Find thalamic channel index and recalculate window length into samples
    ch_index = ch_names.index(thal_chan)
    arousal_duration = int(arousal_duration_seconds * fs)
    
    # Iterate over arousals
    for idx, arousal in arousal_p.iterrows():
        
        # Arousal window
        signal = bp_data[arousal['arousal_start']:arousal['arousal_start']+arousal_duration,ch_index]

        # Signal filtering
        b1, a1 = butter(2, [0.5 / (fs / 2), 4 / (fs / 2)], 'bandpass')
        signal_d = filtfilt(b1, a1, signal)
        
        b2, a2 = butter(3, [4.5 / (fs / 2), 40 / (fs / 2)], 'bandpass')
        signal_tab = filtfilt(b2, a2, signal)
                
        # Power ratio calculation
        result_powerratio = np.median(signal_d ** 2)/np.median(signal_tab ** 2)
                
        # Baseline window
        signal_base = bp_data[arousal['stim_start'] - arousal_duration:arousal['stim_start'],ch_index]

        # Signal filtering
        b1, a1 = butter(2, [0.5 / (fs / 2), 4 / (fs / 2)], 'bandpass')
        signal_d_base = filtfilt(b1, a1, signal_base)
        
        b2, a2 = butter(3, [4.5 / (fs / 2), 40 / (fs / 2)], 'bandpass')
        signal_tab_base = filtfilt(b2, a2, signal_base)
                
        # Power ratio calculation
        result_powerratio_base = np.median(signal_d_base ** 2)/np.median(signal_tab_base ** 2)

        # Results
        power_ratio = pd.concat([power_ratio,pd.DataFrame({
            'patient': patient,
            'time_start': arousal['arousal_start'],
            'channel': thal_chan,
            'power_median': [result_powerratio],
            'power_median_base': [result_powerratio_base]})])

# Relative change computation
power_ratio['relative_change'] = power_ratio['power_median'] / power_ratio['power_median_base']

# Mean value for each arousal
power_ratio_arousals = power_ratio.groupby(['patient','time_start'])['relative_change'].mean().reset_index(drop=False)
    
    
 
    
    
    
    
    
    
    
    
    
