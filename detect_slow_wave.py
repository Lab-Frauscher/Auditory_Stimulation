"""
detect_slow_wave.py

Description:
==================
This script performs automated detection of slow wave events in full-night
European Data Format (EDF) electroencephalographic (EEG) recordings. It extracts
signals using PyEDFlib, dynamically constructs adjacent bipolar montages across
target channels (excluding non-neural auxiliary signals), batches data into 2D matrices,
and executes slow-wave extraction via the SlowWaveDetector pipeline. Detected events
are temporally aligned to absolute UTC Unix timestamps and serialized to disk.

Key Processing Steps:
-------------------
1. Opening EDF recordings via PyEDFlib and resolving global sampling rates and UTC origin timestamps.
2. Standardizing and parsing channel nomenclature across patient cohorts.
3. Constructing bipolar electrode derivation pairs for valid brain channels.
4. Batch reading full-channel signal matrices via PyEDFlib indices.
5. Detecting and filtering slow-wave events (up-state duration: 0.25–1.0 s) using SlowWaveDetector.
6. Aligning event sample positions to absolute UTC Unix timestamps.
7. Exporting structured detection results to serialized Pickle DataFrames.

Dependencies:
-------------------
- pyedflib
- slow_wave_detector
- pandas
- numpy
- importlib
"""

import math
import re
from datetime import timezone
import numpy as np
import pandas as pd
import pyedflib

sw_module = importlib.import_module('slow-wave-detector.slow_wave_detector.slow_wave_detector')
SlowWaveDetector = sw_module.SlowWaveDetector


# --- MAIN PIPELINE EXECUTION ---
patient_ids = [11, 14]

for id_patient in patient_ids:
    file_name = f'demo_data/auditory_stimulation_P{id_patient}_demo.edf'

    # Open EDF recording
    try:
        f = pyedflib.EdfReader(file_name)
    except (OSError, IOError) as e:
        print(f"Error opening EDF file {file_name}: {e}. Skipping patient P{id_patient}.")
        continue

    try:
        unit_conversion = 1.0  # Scale factor for signal amplitude
        fvz = int(f.getSampleFrequency(0))

        # Extract UTC epoch start timestamp with explicit timezone alignment
        start_datetime = f.getStartdatetime()
        if start_datetime is not None:
            time_start = start_datetime.replace(tzinfo=timezone.utc).timestamp()
        else:
            time_start = 0.0

        # --- CHANNEL MAP & NOMENCLATURE STANDARDIZATION ---
        channel_list = [x.strip() for x in f.getSignalLabels()]
        names = [
            i.split('-')[0].split(' ')[1].upper() if id_patient <= 10 else i
            for i in channel_list
        ]
        channel_names = pd.DataFrame({'orig': channel_list, 'name': names})
        channel_names['name'] = channel_names['name'].str.replace(r'^0', 'O', regex=True)

        # Create mapping dictionaries for original labels and PyEDFlib channel indices
        ch_name_to_orig = dict(zip(channel_names['name'], channel_names['orig']))
        orig_to_idx = {label: idx for idx, label in enumerate(channel_list)}

        # --- BIPOLAR MONTAGE CONSTRUCTION ---
        bipolar_montage = []
        norm_names = channel_names['name'].astype(str).tolist()

        for i in range(len(norm_names)):
            curr_name = norm_names[i]
            match = re.match(r"([a-zA-Z]+)([0-9]+)", curr_name)

            if match:
                prefix = match.group(1)
                num = int(match.group(2))

                # Exclude non-brain/auxiliary channels
                if (
                    len(prefix) > 1
                    and prefix.upper() not in ["EKG", "ECG", "DC", "FP"]
                ):
                    target_num = num + 1
                    for j in range(len(norm_names)):
                        next_name = norm_names[j]
                        next_match = re.match(r"([a-zA-Z]+)([0-9]+)", next_name)

                        if next_match:
                            if (
                                next_match.group(1) == prefix
                                and int(next_match.group(2)) == target_num
                            ):
                                bipolar_montage.append(f"{curr_name}-{next_name}")
                                break

        results = pd.DataFrame([])

        # --- BATCHED SIGNAL EXTRACTION & SLOW WAVE DETECTION ---
        chunk_size = 10
        num_chunks = math.ceil(len(bipolar_montage) / chunk_size)

        for curr_channels in np.array_split(bipolar_montage, num_chunks):
            if len(curr_channels) == 0:
                continue

            signals = []
            valid_channels = []

            for name in curr_channels:
                name1, name2 = name.split('-')[0], name.split('-')[1]
                if name1 in ch_name_to_orig and name2 in ch_name_to_orig:
                    orig_name1 = ch_name_to_orig[name1]
                    orig_name2 = ch_name_to_orig[name2]

                    idx1 = orig_to_idx[orig_name1]
                    idx2 = orig_to_idx[orig_name2]

                    # Read channel signals directly via PyEDFlib indices
                    data1 = f.readSignal(idx1)
                    data2 = f.readSignal(idx2)

                    # Derive bipolar differential signal
                    signal = -1.0 * (data1 - data2) * unit_conversion
                    signals.append(signal)
                    valid_channels.append(name)

            if not signals:
                print(f"Warning: No valid channel pairs found in current batch for patient P{id_patient}.")
                continue

            # Construct 2D array (N_channels x N_samples) and handle NaNs
            data_2d = np.array(signals)
            data_2d = np.nan_to_num(data_2d, nan=0.0)

            # Run Slow Wave Detector algorithm
            detector = SlowWaveDetector(eeg_data=data_2d, fs=fvz)
            sw_detections_all = detector.extract_slow_wave_events()
            sw_detections = detector.filter_events(
                sw_detections_all,
                min_up_state_duration=0.25,
                max_up_state_duration=1.0
            )

            # --- PROCESS DETECTED EVENTS ---
            for chan_result in range(len(sw_detections)):
                result_dict = [event.to_dict() for event in sw_detections[chan_result]]

                if len(result_dict) > 0:
                    curr_result = pd.DataFrame(result_dict)
                    curr_result['unix_start'] = (
                        curr_result['start_sample_1_based'] / curr_result['sampling_frequency']
                    ) + time_start
                    curr_result['unix_end'] = (
                        curr_result['end_sample_1_based'] / curr_result['sampling_frequency']
                    ) + time_start
                    curr_result['channel_name'] = curr_result['channel_idx'].apply(
                        lambda x: valid_channels[x - 1]
                    )
                    curr_result = curr_result[['channel_name', 'unix_start', 'unix_end']]

                    results = pd.concat([results, curr_result], ignore_index=True)

        # Save cumulative results per patient
        #results.to_pickle(f'detector_results/p{id_patient}_slow-wave.pkl')

    finally:
        f.close()
