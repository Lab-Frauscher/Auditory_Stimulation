"""
hfo_ied_detection_seeg.py

Description:
==================
This script performs automated, parallelized detection of High-Frequency Oscillations
(HFOs) and Interictal Epileptiform Discharges (IEDs / Spikes) in stereo-EEG (SEEG)
recordings. Utilizing European Data Format (EDF) input files via PyEDFlib, it
dynamically constructs consecutive bipolar derivation montages, resamples signals to
a uniform sampling rate of 2000 Hz, and evaluates 60-second non-overlapping signal
windows using algorithmically validated detectors (Nicolas HFO detector and Janca spike
detector). Detected events undergo spatiotemporal artifact filtering to remove
multichannel volume-conducted artifacts and intra-channel temporal clustering before
exporting absolute Unix timestamps.

Key Processing Steps:
-------------------
1. Reading EDF metadata and standardizing channel nomenclature across patient cohorts.
2. Dynamic construction of adjacent bipolar channel pairs (excluding non-neural EKG/ECG/DC/FP signals).
3. Multiprocessing distribution of channel pairs across CPU workers.
4. Signal extraction, bipolar differential derivation, and polyphase resampling to 2000 Hz.
5. Windowed detection of HFOs (detect_hfo_nicolas) and IEDs (detect_spikes_janca).
6. Remapping sample-level event detections back to the native EDF sampling frequency.
7. Post-processing artifact rejection:
   - Spatial filtering: Discarding widespread events (> 50% channels within a 50 ms window).
   - Temporal filtering: Enforcing a 300 ms refractory period between consecutive spikes per channel.
8. Serialization of validated HFO and IED event DataFrames with absolute UTC Unix timestamps.

Dependencies:
-------------------
- pyedflib
- epycom (epycom.event_detection.hfo, epycom.event_detection.spike)
- scipy
- pandas
- numpy
"""

import math
import re
from multiprocessing import Pool
from datetime import timezone
import numpy as np
import pandas as pd
import pyedflib
from epycom.event_detection.hfo.nicolas_detector import detect_hfo_nicolas
from epycom.event_detection.spike.janca_detector import detect_spikes_janca
from scipy import signal


def process_channel_hfo_and_ied(channel_name):
    """
    Extracts bipolar SEEG signals, resamples data, and computes HFO and IED detections.

    Parameters:
    -----------
    channel_name : list of int
        Tuple/list containing indices of the primary and reference channels [ch1_idx, ch2_idx].

    Returns:
    --------
    hfo_channel_df : pandas.DataFrame
        Detected HFO events containing sample bounds, oscillation frequency, and metadata.
    ied_channel_df : pandas.DataFrame
        Detected IED events containing sample peak indices and channel identifiers.
    """
    # --- 1. SIGNAL EXTRACTION AND BIPOLAR DERIVATION ---
    # Extract raw signal vectors for specified channel index pair
    data1 = f.readSignal(channel_name[0])
    data2 = f.readSignal(channel_name[1])

    data = -1 * np.array([data1, data2])
    data_bipolar = data[0] - data[1]

    # Resample signal to standard 2000 Hz using high-efficiency polyphase filtering if needed
    if fvz != 2000:
        data_bipolar = signal.resample_poly(data_bipolar, 2000, fvz)
        fvz_channel = 2000
    else:
        fvz_channel = 2000

    # Format channel label for output structure
    ch_label = f"{channel_names['name'][channel_name[0]]}-{channel_names['name'][channel_name[1]]}"

    # Define non-overlapping 60-second window length (samples)
    window_size = 60 * fvz_channel

    # --- 2. HIGH-FREQUENCY OSCILLATION (HFO) DETECTION ---
    hfo_list = []
    for start in range(0, len(data_bipolar) - window_size + 1, window_size):
        window = data_bipolar[start: start + window_size]
        res_hfo = detect_hfo_nicolas(window, fs=fvz_channel)

        res_df = pd.DataFrame(
            res_hfo, columns=['event_start', 'event_stop', 'osc_frequency']
        )
        # Remap sample indices back to native EDF sampling frequency scale
        res_df['event_start'] = (
                (start + res_df['event_start']) * (fvz / fvz_channel)
        ).astype(int)
        res_df['event_stop'] = (
                (start + res_df['event_stop']) * (fvz / fvz_channel)
        ).astype(int)
        res_df.insert(0, 'ID', id_patient)
        res_df.insert(1, 'channel', ch_label)
        hfo_list.append(res_df)

    hfo_list = [df for df in hfo_list if not df.empty]
    hfo_channel_df = (
        pd.concat(hfo_list, ignore_index=True) if hfo_list else pd.DataFrame()
    )

    # --- 3. INTERICTAL EPILEPTIFORM DISCHARGE (IED) DETECTION ---
    ied_list = []
    for start in range(0, len(data_bipolar) - window_size + 1, window_size):
        window = data_bipolar[start: start + window_size]
        res_ied = detect_spikes_janca(window, fs=fvz_channel)

        ied_df = pd.DataFrame({
            'ID': id_patient,
            'channel': [ch_label] * len(res_ied),
            'sample_peak': [
                int((start + row[0]) * (fvz / fvz_channel)) for row in res_ied
            ],
        })
        ied_list.append(ied_df)

    ied_list = [df for df in ied_list if not df.empty]
    ied_channel_df = (
        pd.concat(ied_list, ignore_index=True) if ied_list else pd.DataFrame()
    )
    ied_channel_df = ied_channel_df.drop_duplicates()

    return hfo_channel_df, ied_channel_df


def ied_postprocessing(
        ied_detection, channel_number, window_cross_s=0.05, window_same_s=0.3
):
    """
    Applies spatiotemporal filtering to reject multichannel artifacts and temporal clusters.

    Parameters:
    -----------
    ied_detection : pandas.DataFrame
        Raw spike detections containing timestamps and channel identifiers.
    channel_number : int
        Total number of evaluated bipolar channels.
    window_cross_s : float, default=0.05
        Coincidence window (seconds) for detecting widespread volume-conducted artifacts.
    window_same_s : float, default=0.3
        Refractory period (seconds) for rejecting consecutive spikes in the same channel.

    Returns:
    --------
    pandas.DataFrame
        Filtered IED detection dataset.
    """
    if ied_detection.empty:
        return ied_detection

    ied_clean = (
        ied_detection.copy().sort_values('time_peak').reset_index(drop=True)
    )
    times = ied_clean['time_peak'].values
    channels = ied_clean['channel'].values

    # --- STEP 1: MULTICHANNEL ARTIFACT REJECTION ---
    # Identify coincident spikes across > 50% of channels within the cross-window (50 ms)
    left_idxs = np.searchsorted(times, times - window_cross_s / 2, side='left')
    right_idxs = np.searchsorted(
        times, times + window_cross_s / 2, side='right'
    )

    to_drop_artifact = np.zeros(len(ied_clean), dtype=bool)
    for i in range(len(ied_clean)):
        chans_in_window = np.unique(channels[left_idxs[i]: right_idxs[i]])
        if len(chans_in_window) > (channel_number / 2):
            to_drop_artifact[i] = True

    ied_clean = (
        ied_clean[~to_drop_artifact]
        .sort_values('time_peak')
        .reset_index(drop=True)
    )

    # --- STEP 2: INTRA-CHANNEL REFRACTORY FILTERING ---
    # Enforce minimum temporal spacing (300 ms) between consecutive spikes per channel
    keep_mask = np.zeros(len(ied_clean), dtype=bool)
    for chan in ied_clean['channel'].unique():
        chan_indices = ied_clean.index[ied_clean['channel'] == chan].tolist()
        if not chan_indices:
            continue

        chan_times = ied_clean.loc[chan_indices, 'time_peak'].values
        keep_mask[chan_indices[0]] = True
        last_kept_time = chan_times[0]

        for i in range(1, len(chan_times)):
            current_time = chan_times[i]
            if (current_time - last_kept_time) >= window_same_s:
                keep_mask[chan_indices[i]] = True
                last_kept_time = current_time

    return ied_clean[keep_mask].sort_values('time_peak').reset_index(drop=True)


# --- MAIN PIPELINE EXECUTION ---
patients = [11, 14]

for id_patient in patients:
    file_name = f'demo_data/auditory_stimulation_P{id_patient}_demo.edf'

    # Open recording via PyEDFlib
    f = pyedflib.EdfReader(file_name)

    # Extract global sampling frequency and epoch origin timestamp (UTC aligned)
    fvz = int(f.getSampleFrequency(0))

    start_datetime = f.getStartdatetime()
    if start_datetime is not None:
        time_start = start_datetime.replace(tzinfo=timezone.utc).timestamp()
    else:
        time_start = 0.0

    # Parse and normalize channel nomenclature across recording formats
    channel_list = f.getSignalLabels()
    names = [
        i.split('-')[0].split(' ')[1].upper() if id_patient <= 10 else i
        for i in channel_list
    ]
    channel_names = pd.DataFrame({'orig': channel_list, 'name': names})
    channel_names['name'] = channel_names['name'].str.replace(
        r'^0', 'O', regex=True
    )

    # --- CONSTRUCT BIPOLAR DERIVATION MONTAGE ---
    # Pair adjacent contact numbers within electrode shafts, excluding non-neural channels
    bipolar_montage_indices = []
    indices = channel_names.index.tolist()
    names_list = channel_names['name'].astype(str).tolist()

    for i in range(len(indices)):
        curr_name = names_list[i]
        match = re.match(r'([a-zA-Z]+)([0-9]+)', curr_name)

        if match:
            prefix = match.group(1)
            num = int(match.group(2))

            # Exclude auxiliary artifact channels (EKG, ECG, DC, FP)
            if (
                    len(prefix) > 1
                    and prefix.upper() not in ['EKG', 'ECG', 'DC', 'FP']
            ):
                target_num = num + 1
                for j in range(len(indices)):
                    next_match = re.match(
                        r'([a-zA-Z]+)([0-9]+)', names_list[j]
                    )
                    if (
                            next_match
                            and next_match.group(1) == prefix
                            and int(next_match.group(2)) == target_num
                    ):
                        bipolar_montage_indices.append(
                            [indices[i], indices[j]]
                        )
                        break

    # --- PARALLEL PROCESSING VIA PROCESS POOL ---
    N_processes = 15

    all_hfo_results = []
    all_ied_results = []

    # Partition channel pairs into batches for parallel pool mapping
    chunks = np.array_split(
        bipolar_montage_indices,
        math.ceil(len(bipolar_montage_indices) / N_processes),
    )

    for iter_map_i in chunks:
        if N_processes > 1:
            with Pool(N_processes) as mp:
                results = mp.map(process_channel_hfo_and_ied, iter_map_i)
        else:
            results = [
                process_channel_hfo_and_ied(pair) for pair in iter_map_i
            ]

        for hfo_res, ied_res in results:
            all_hfo_results.append(hfo_res)
            all_ied_results.append(ied_res)

    # --- EXPORT HFO METRICS ---
    hfo_detection = (
        pd.concat(all_hfo_results, ignore_index=True)
        if all_hfo_results
        else pd.DataFrame()
    )
    if not hfo_detection.empty:
        hfo_detection['unix_start'] = time_start + (
                hfo_detection['event_start'].astype('float64') / fvz
        )
        hfo_detection['unix_stop'] = time_start + (
                hfo_detection['event_stop'].astype('float64') / fvz
        )
        hfo_detection.to_pickle(f'detector_results/p{id_patient}_hfo.pkl')

    # --- POST-PROCESS AND EXPORT IED METRICS ---
    ied_detection = (
        pd.concat(all_ied_results, ignore_index=True)
        if all_ied_results
        else pd.DataFrame()
    )
    if not ied_detection.empty:
        ied_detection['time_peak'] = time_start + (
                ied_detection['sample_peak'].astype('float64') / fvz
        )

        # Apply spatiotemporal artifact cleaning
        clear_ied = ied_postprocessing(
            ied_detection, len(bipolar_montage_indices)
        )
        clear_ied['ID'] = f'P{id_patient}'

        clear_ied.to_pickle(f'detector_results/p{id_patient}_ied.pkl')

    f.close()
    print(f'Patient P{id_patient} processing complete.')
