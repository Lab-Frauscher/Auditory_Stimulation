"""
ied_slow-wave_power-ratio_analysis.py

Description:
==================
This script implements a feature extraction pipeline to evaluate interictal 
epileptiform discharge (IED), slow-wave co-occurrences, and scalp EEG power 
ratios across individual stereo-EEG (SEEG) channels during sleep arousals. 
The analysis measures raw spike counts within paired 3-second baseline and 
post-event epochs across experimental conditions (evoked, spontaneous, and control). 
Each channel-level observation is integrated with slow-wave detection status, pre-stimulus 
scalp power ratio metrics, and clinical/anatomical metadata (SOZ vs. non-SOZ, 
mesiotemporal vs. neocortical).

Key Processing Steps:
-------------------
1. Loading arousal event annotations, SEEG contact metadata, active channel flags, 
   and pre-computed scalp power ratios.
2. Constructing patient-specific bipolar channel maps and mapping metadata attributes.
3. Loading automated IED spike timestamps and slow-wave event boundaries.
4. Defining condition-dependent time windows for pre-stimulus baseline (-3 to 0 s), 
   post-arousal event (+0 to +3 s), and slow-wave evaluation windows.
5. Evaluating slow-wave co-occurrence and counting unclustered raw IED spikes per channel.
6. Extracting corresponding time-aligned pre-stimulus scalp power ratio values.
7. Aggregating channel-level results and executing natural sorting across patients and events.

Dependencies:
-------------------
- numpy
- pandas
- openpyxl
"""

import pandas as pd
import numpy as np

# --- 1. DATA LOADING AND INITIALIZATION ---
# Load arousal event annotations and create global event index
arousals = pd.read_excel('arousal_annotations.xlsx')
arousals['Event number'] = range(1, len(arousals) + 1)

# Load electrode metadata and normalize contact labels
channel_labels = pd.read_excel('channel_labels.xlsx')
channel_labels['Contact'] = channel_labels['Contact'].str.upper()

# Load scalp power ratio data and linearly interpolate missing values per patient
power_ratio_df = pd.read_pickle('powerratio_scalp_stimulation.pkl')
power_ratio_df['power_ratio'] = power_ratio_df.groupby('ID')['power_ratio'].transform(
    lambda x: x.interpolate(method='linear', limit_direction='both')
)

# Load pre-computed active channel flags
active_channels_df = pd.read_pickle('active_channels.pkl')

# Construct fast lookup map for channel activity status: (Patient_ID, channel) -> 1/0
active_map = {}
for _, row in active_channels_df.iterrows():
    active_map[(row['ID'], row['channel'])] = row['Active channel']

# --- 2. HELPER FUNCTIONS ---
def get_anatomical_and_soz(bipolar_chan, loc_dict, soz_dict):
    """
    Determines concordant anatomical region and clinical zone classifications for a bipolar channel.

    Parameters:
    -----------
    bipolar_chan : str
        Bipolar channel identifier (e.g., 'A1-A2').
    loc_dict : dict
        Mapping of contact labels to anatomical locations ('mesiotemporal' vs 'neocortical').
    soz_dict : dict
        Mapping of contact labels to clinical zones ('SOZ' vs 'non-SOZ').

    Returns:
    --------
    loc_group : str
        Anatomical classification ('mesiotemporal', 'neocortical', or 'mixed/unknown').
    soz_group : str
        Clinical zone classification ('SOZ', 'non-SOZ', or 'mixed/unknown').
    """
    c1, c2 = [p.strip() for p in str(bipolar_chan).split('-')]
    loc1, loc2 = str(loc_dict.get(c1, '')).lower(), str(loc_dict.get(c2, '')).lower()
    soz1, soz2 = str(soz_dict.get(c1, '')).lower(), str(soz_dict.get(c2, '')).lower()
    
    loc_group = 'mesiotemporal' if loc1 == 'mesiotemporal' and loc2 == 'mesiotemporal' else ('neocortical' if loc1 == 'neocortical' and loc2 == 'neocortical' else 'mixed/unknown')
    soz_group = 'SOZ' if soz1 == 'soz' and soz2 == 'soz' else ('non-SOZ' if soz1 == 'non-soz' and soz2 == 'non-soz' else 'mixed/unknown')
    
    return loc_group, soz_group

# --- 3. PATIENT AND EVENT PROCESSING PIPELINE ---
all_results = []

for patient_id in ['P11', 'P14']:
    # Load patient-specific automated IED and slow-wave detection results
    IED_id = pd.read_pickle('detector_results/' + patient_id.lower() + '_ied.pkl')
    slow_waves_id = pd.read_pickle('detector_results/' + patient_id.lower() + '_slow-wave.pkl')
    
    # Filter and chronologically sort arousal events for the current patient
    arousal_id = arousals.loc[arousals['ID'] == patient_id].sort_values(by=['utc_stim_start', 'utc_arousal_start'])
    
    # Restrict analysis to events within a 4-hour post-stimulation protocol window
    fragmentation_stop = arousal_id['utc_stim_start'].max()
    arousal_id = arousal_id.loc[(arousal_id['utc_arousal_start'] <= fragmentation_stop + 4 * 3600) | (arousal_id['utc_arousal_start'].isna())]
    
    # Filter channel labels for current patient
    channel_label_id = channel_labels.loc[channel_labels['Patient ID'] == patient_id]
    
    # Extract channel metadata lookup dictionaries
    loc_dict = dict(zip(channel_label_id['Contact'].astype(str).str.upper().str.strip(), channel_label_id['mesiotemporal vs neocortical'].astype(str).str.strip()))
    soz_dict = dict(zip(channel_label_id['Contact'].astype(str).str.upper().str.strip(), channel_label_id['SOZ vs non-SOZ'].astype(str).str.strip()))
    
    # Reconstruct all adjacent bipolar contact pairs belonging to the same electrode stem
    contacts = channel_label_id['Contact'].astype(str).tolist()
    unique_bipolar_chans = []
    for i in range(len(contacts) - 1):
        prefix1 = ''.join([c for c in contacts[i] if c.isalpha()])
        prefix2 = ''.join([c for c in contacts[i+1] if c.isalpha()])
        if prefix1 == prefix2 and prefix1 != '':
            unique_bipolar_chans.append(f"{contacts[i]}-{contacts[i+1]}")
    
    patient_spikes_df = IED_id[['time_peak', 'channel']]
    
    # --- AROUSAL EVENT ITERATION ---
    for index, row in arousal_id.iterrows():
        a_type = row['type']
        stim_time = row['utc_stim_start']
        arousal_start = row['utc_arousal_start']
        ev_num = row['Event number']
        
        # Define 3-second baseline (pre), event (post), and slow-wave search windows
        if a_type == 'evoked':
            pre_start, pre_end = stim_time - 3, stim_time
            post_start, post_end = arousal_start, arousal_start + 3
            sw_win_start = arousal_start if not pd.isna(arousal_start) else 0
            sw_win_end = sw_win_start + row.get('duration', 0)
            
        elif a_type == 'control':
            pre_start, pre_end = stim_time - 3, stim_time
            post_start, post_end = stim_time + 1.076, stim_time + 1.076 + 3
            sw_win_start = stim_time if not pd.isna(stim_time) else 0
            sw_win_end = sw_win_start + 6.5
            
        elif a_type == 'spontaneous':
            pre_start, pre_end = arousal_start - 3, arousal_start
            post_start, post_end = arousal_start, arousal_start + 3
            sw_win_start = arousal_start if not pd.isna(arousal_start) else 0
            sw_win_end = sw_win_start + row.get('duration', 0)
        else:
            continue
            
        # Determine slow-wave co-occurrence within the defined window (binary flag: 1=Yes, 0=No)
        sw_present = int(any((slow_waves_id['unix_start'] <= sw_win_end) & (slow_waves_id['unix_end'] >= sw_win_start)))
        
        # Filter raw spike peak timestamps within baseline and post-event windows
        pre_spikes = patient_spikes_df[(patient_spikes_df['time_peak'] >= pre_start) & (patient_spikes_df['time_peak'] <= pre_end)]
        post_spikes = patient_spikes_df[(patient_spikes_df['time_peak'] >= post_start) & (patient_spikes_df['time_peak'] <= post_end)]
        
        # Compute unclustered raw spike counts per channel
        pre_counts = pre_spikes['channel'].value_counts().to_dict()
        post_counts = post_spikes['channel'].value_counts().to_dict()
        
        # Retrieve time-aligned pre-stimulus scalp power ratio
        row_pr = power_ratio_df.loc[np.isclose(power_ratio_df['utc_stim_start'], stim_time, atol=1e-4)]
        
        # --- CHANNEL-LEVEL DATA AGGREGATION ---
        for chan in unique_bipolar_chans:
            loc_group, soz_group = get_anatomical_and_soz(chan, loc_dict, soz_dict)
            is_active = active_map.get((patient_id, chan), 0)
            
            pr_val = row_pr['power_ratio'].values[0] if (not row_pr.empty and 'power_ratio' in row_pr.columns) else np.nan
            
            result_row = {
                'Patient ID': patient_id,
                'Event type': a_type,
                'Event number': ev_num,
                'Stimulus name': row.get('stim_intensity', 'N/A'),
                'Pre-stimulus Sleep Stage': row.get('stage', np.nan),
                'Post-stimulus Sleep Stage': row.get('stage_post', np.nan),
                'SEEG channel': chan,
                'Active channel (1=Y, 0=N)': int(is_active),
                'Slow-wave presence (1=Y, 0=N)': sw_present,
                'Pre-stimulus spectral power ratio (scalp)': pr_val,
                'Pre-stimulus IED no. on contact': pre_counts.get(chan, 0),
                'Arousal (/event) IED no. on contact': post_counts.get(chan, 0),
                'Channel type (mesiotemporal or neocortical)': loc_group,
                'Channel type (SOZ or non-SOZ)': soz_group
            }
            all_results.append(result_row)

# --- 4. DATA SYNTHESIS AND NATURAL SORTING ---
final_df = pd.DataFrame(all_results)

# Apply natural sorting across Patient IDs, Event Types, Event Numbers, and Channel Labels
final_df = final_df.sort_values(
    ['Patient ID', 'Event type', 'Event number', 'SEEG channel'],
    key=lambda col: col.str.replace('P', '').astype(int) if col.name == 'Patient ID' else col
)