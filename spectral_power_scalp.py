"""
spectral_power_scalp.py

Description:
==================
This script extracts pre-stimulus scalp EEG epochs from European Data Format (EDF) 
recordings to calculate spectral power ratios across auditory stimulation events.
It dynamically selects available bipolar scalp channels (prioritizing FZ-CZ, F3-C3, 
or CZ-PZ), segments 3-second baseline windows immediately prior to stimulation onsets, 
applies bandpass Butterworth filters, and computes the ratio of delta-band power (0.5-4 Hz) 
to higher-frequency power (4.5-30 Hz).

Dependencies:
-------------------
- pyedflib
- scipy
- pandas
- numpy
"""

import numpy as np
import pandas as pd
import pyedflib
from datetime import timezone
from scipy.signal import butter, sosfiltfilt

# --- 1. DATA INITIALIZATION AND ANNOTATIONS LOADING ---
annot = pd.read_excel('arousal_annotations.xlsx')

powerratio_results = pd.DataFrame({
    'ID': [],
    'utc_stim_start': [],
    'power_ratio': [],
    'channel': [],
    'type': []
})

target_patients = [11, 14]

# --- 2. PATIENT RECORDING PROCESSING PIPELINE ---
for id_pacient in target_patients:
    file_name = f'demo_data/auditory_stimulation_P{id_pacient}_demo.edf'

    # Load raw EDF recording using pyedflib
    f = pyedflib.EdfReader(file_name)

    try:
        unit_conversion = 1.0  # Unit scale factor
        
        # Calculate UTC epoch origin timestamp for time-alignment
        start_datetime = f.getStartdatetime()
        if start_datetime is not None:
            time_start = start_datetime.replace(tzinfo=timezone.utc).timestamp()
        else:
            time_start = 0.0

        # --- 3. CHANNEL SELECTION AND FALLBACK MAPPING ---
        channel_list = f.getSignalLabels()
        names = [i.split('-')[0].split(' ')[1].upper() if id_pacient <= 10 else i for i in channel_list]
        channel_names = pd.DataFrame({
            'orig': channel_list, 
            'name': names,
            'idx': range(len(channel_list))
        })
        channel_names['name'] = channel_names['name'].str.replace(r'^0', 'O', regex=True)

        # Priority hierarchy for selecting bipolar scalp channels
        if channel_names['name'].str.contains('FZ').any() and channel_names['name'].str.contains('CZ').any():
            channel = 'FZ-CZ'
            ch0_idx = channel_names.loc[channel_names['name'] == 'FZ']['idx'].item()
            ch1_idx = channel_names.loc[channel_names['name'] == 'CZ']['idx'].item()
        elif channel_names['name'].str.contains('F3').any() and channel_names['name'].str.contains('C3').any():
            channel = 'F3-C3'
            ch0_idx = channel_names.loc[channel_names['name'] == 'F3']['idx'].item()
            ch1_idx = channel_names.loc[channel_names['name'] == 'C3']['idx'].item()
        elif channel_names['name'].str.contains('CZ').any() and channel_names['name'].str.contains('PZ').any():
            channel = 'CZ-PZ'
            ch0_idx = channel_names.loc[channel_names['name'] == 'CZ']['idx'].item()
            ch1_idx = channel_names.loc[channel_names['name'] == 'PZ']['idx'].item()
        else:
            print(f"P{id_pacient}: Missing required scalp channel pairs. Skipping patient.")
            continue

        # Sampling rate (Hz) derived from selected channel
        fvz = int(f.getSampleFrequency(ch0_idx))

        # Filter events for current patient with valid stimulation start timestamps
        stim_p = annot.loc[(annot['ID'] == f'P{id_pacient}') & ~(annot['utc_stim_start'].isna())]

        # --- 4. PRE-STIMULUS EPOCH FILTERING & FEATURE EXTRACTION ---
        for idx, stim in stim_p.iterrows():
            # Convert absolute UTC timestamp to relative recording duration (seconds)
            tmax = stim['utc_stim_start'] - time_start
            tmin = tmax - 3  # 3-second baseline window prior to stimulus onset

            # Convert time window (s) to sample indices for pyedflib
            start_sample = int(np.round(tmin * fvz))
            n_samples = int(np.round((tmax - tmin) * fvz))

            # Read signal slices using index and sample length
            data_ch0 = f.readSignal(ch0_idx, start_sample, n_samples)
            data_ch1 = f.readSignal(ch1_idx, start_sample, n_samples)

            data = np.array([data_ch0, data_ch1])
            data = -1 * data * unit_conversion

            # Construct single bipolar differential channel signal
            signal = data[0] - data[1]

            # Filter 1: Delta frequency band (0.5 - 4.0 Hz, 2nd order Butterworth, zero-phase)
            sos_delta = butter(2, [0.5, 4.0], 'bandpass', fs=fvz, output='sos')
            signal_d = sosfiltfilt(sos_delta, signal)

            # Filter 2: Theta/Alpha/Beta frequency band (4.5 - 30.0 Hz, 3rd order Butterworth, zero-phase)
            sos_tab = butter(3, [4.5, 30.0], 'bandpass', fs=fvz, output='sos')
            signal_tab = sosfiltfilt(sos_tab, signal)

            # Feature Computation: Ratio of median squared power (Delta / Higher Bands)
            result_powerratio = np.median(signal_d ** 2) / np.median(signal_tab ** 2)

            # Aggregate row entry
            event_df = pd.DataFrame({
                'ID': [f'P{id_pacient}'],
                'utc_stim_start': [stim['utc_stim_start']],
                'power_ratio': [result_powerratio],
                'channel': [channel],
                'type': [stim['type']]
            })
            powerratio_results = pd.concat([powerratio_results, event_df], ignore_index=True)

    finally:
        f.close()

# --- 5. EXPORT / STORAGE (OPTIONAL) ---
# powerratio_results.to_pickle('powerratio_scalp_stimulation.pkl')
