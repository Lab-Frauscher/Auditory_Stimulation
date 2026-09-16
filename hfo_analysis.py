"""
hfo_analysis.py

Description:
==================
This script implements a signal feature extraction and event-merging pipeline to 
quantify High-Frequency Oscillation (HFO) rates across individual intracranial EEG (SEEG) 
channels during sleep arousals. The analysis evaluates HFO occurrences (>= 80 Hz) within 
paired 3-second pre- and post-event windows across different experimental conditions 
(evoked, spontaneous, and control). Overlapping HFO detections are temporally merged 
using a reciprocal overlap threshold (>= 50%). Results are formatted per channel and 
annotated with anatomical (mesiotemporal vs. neocortical) and clinical (SOZ vs. non-SOZ) 
classifications.

Key Processing Steps:
-------------------
1. Loading arousal event annotations, contact metadata, and active channel lookup tables.
2. Generating unique global event identifiers and constructing patient-specific bipolar channel maps.
3. Loading patient HFO detection records and filtering for ripple/fast-ripple frequency bands (>= 80 Hz).
4. Window segmentation for baseline/pre-stimulus (-3 to 0 s) and post-arousal/post-stimulus (+0 to +3 s) conditions.
5. Merging temporally overlapping HFO detections (>= 50% overlap relative to shorter event) per channel.
6. Mapping channel-level HFO counts alongside anatomical, clinical, and signal activity metadata.
7. Aggregating structured outputs with natural sorting across patient IDs and event numbers.

Dependencies:
-------------------
- numpy
- pandas
- openpyxl
"""

import os
import pandas as pd
import numpy as np

# --- 1. DATA LOADING AND INITIALIZATION ---
# Load arousal event annotations and SEEG channel metadata
arousals = pd.read_excel('arousal_annotations.xlsx')
channel_labels = pd.read_excel('channel_labels.xlsx')
channel_labels['Contact'] = channel_labels['Contact'].str.upper()

# Assign unique global identifiers across all experimental events
arousals['Event number'] = range(1, len(arousals) + 1)

# Load pre-computed active channel flags
active_channels_df = pd.read_pickle('active_channels.pkl')

# Construct fast lookup map for channel activity status: (Patient_ID, channel) -> 1/0
active_map = {}
for _, row in active_channels_df.iterrows():
    active_map[(row['ID'], row['channel'])] = row['Active channel']

all_results = []
unique_patients = arousals['ID'].unique()

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

def count_merged_hfos(df_window):
    """
    Iterates through detected HFOs within a time window, sorts them chronologically, 
    and merges events exhibiting >= 50% temporal overlap relative to the shorter duration.

    Parameters:
    -----------
    df_window : pd.DataFrame
        DataFrame containing detected HFO events with 'channel', 'unix_start', and 'unix_stop' fields.

    Returns:
    --------
    pd.Series
        Series containing counts of unique merged HFO events per channel.
    """
    if df_window.empty:
        return pd.Series(dtype=int)
        
    counts = {}
    for channel, group in df_window.groupby('channel'):
        # Sort detected events chronologically by onset time
        events = group.sort_values('unix_start').to_dict('records')
        
        if not events:
            counts[channel] = 0
            continue
            
        merged_count = 0
        current_event = events[0]
        
        for next_event in events[1:]:
            s1 = current_event['unix_start']
            e1 = current_event['unix_stop'] 
            
            s2 = next_event['unix_start']
            e2 = next_event['unix_stop']
            
            # Calculate temporal overlap and event durations
            overlap = max(0, min(e1, e2) - max(s1, s2))
            len1 = e1 - s1
            len2 = e2 - s2
            min_len = min(len1, len2)
            
            # Merge events if overlap exceeds 50% of the shorter event duration
            if min_len > 0 and (overlap / min_len) >= 0.5:
                current_event['unix_stop'] = max(e1, e2)
            else:
                merged_count += 1
                current_event = next_event
                
        # Include final event in current window
        merged_count += 1 
        counts[channel] = merged_count
        
    return pd.Series(counts)

# --- 3. PATIENT AND EVENT PROCESSING PIPELINE ---
for patient_id in unique_patients:
    p_num = patient_id.lower().replace('p', '')
    hfo_path = f'detector_results/p{p_num}_hfo.pkl'
    
    if not os.path.exists(hfo_path): 
        continue
    hfo_data = pd.read_pickle(hfo_path)
    
    # Filter detections to target HFO frequency band (>= 80 Hz)
    hfo_data = hfo_data[hfo_data['osc_frequency'] >= 80]
    
    # Extract channel metadata mappings for current patient
    p_labels = channel_labels[channel_labels['Patient ID'] == patient_id]
    
    loc_dict = dict(zip(p_labels['Contact'].astype(str).str.upper().str.strip(), p_labels['mesiotemporal vs neocortical'].astype(str).str.strip()))
    soz_dict = dict(zip(p_labels['Contact'].astype(str).str.upper().str.strip(), p_labels['SOZ vs non-SOZ'].astype(str).str.strip()))
    
    # Reconstruct all adjacent bipolar contact pairs for current patient
    contacts = p_labels['Contact'].astype(str).tolist()
    unique_chans = []
    for i in range(len(contacts) - 1):
        prefix1 = ''.join([c for c in contacts[i] if c.isalpha()])
        prefix2 = ''.join([c for c in contacts[i+1] if c.isalpha()])
        if prefix1 == prefix2 and prefix1 != '':
            unique_chans.append(f"{contacts[i]}-{contacts[i+1]}")
    
    # Apply protocol time limit on arousal events (within 4 hours post-fragmentation)
    p_arousals = arousals[arousals['ID'] == patient_id].copy()
    fragmentation_stop = p_arousals['utc_stim_start'].max()
    p_arousals = p_arousals.loc[
        (p_arousals['utc_arousal_start'] <= fragmentation_stop + 4 * 3600) | 
        (p_arousals['utc_arousal_start'].isna())
    ]
    
    # Process individual arousal events
    for _, row in p_arousals.iterrows():
        a_type = row['type']
        t_stim = row['utc_stim_start']
        t_arou = row['utc_arousal_start']
        ev_num = row['Event number']
        stim_intensity = row.get('stim_intensity', np.nan)
        
        sleep_stage_pre = row.get('stage', np.nan)
        sleep_stage_post = row.get('stage_post', np.nan)
        
        # Define baseline (pre) and event (post) analysis windows based on event condition
        if a_type == 'evoked':
            w_pre, w_post = (t_stim - 3, t_stim), (t_arou, t_arou + 3)
        elif a_type == 'control':
            w_pre, w_post = (t_stim - 3, t_stim), (t_stim + 1.076, t_stim + 1.076 + 3)
        elif a_type == 'spontaneous':
            w_pre, w_post = (t_arou - 3, t_arou), (t_arou, t_arou + 3)
        else: 
            continue

        # Extract HFO detections within defined time windows
        hfo_data = hfo_data.loc[hfo_data['osc_frequency'] >= 80]
        hfo_pre = hfo_data[(hfo_data['unix_start'] >= w_pre[0]) & (hfo_data['unix_start'] <= w_pre[1])]
        hfo_post = hfo_data[(hfo_data['unix_start'] >= w_post[0]) & (hfo_data['unix_start'] <= w_post[1])]
        
        # Merge overlapping detections and compute channel event counts
        hfo_pre_counts = count_merged_hfos(hfo_pre)
        hfo_post_counts = count_merged_hfos(hfo_post)
        
        # Compile channel-level results
        for chan in unique_chans:
            loc_group, soz_group = get_anatomical_and_soz(chan, loc_dict, soz_dict)
            
            is_active = active_map.get((patient_id, chan), 0)
            
            all_results.append({
                'Patient ID': patient_id,
                'Event type': a_type,
                'Event number': ev_num,
                'Stimulus name': stim_intensity,
                'Pre-stimulus Sleep Stage': sleep_stage_pre,
                'Post-stimulus Sleep Stage': sleep_stage_post,
                'SEEG channel': chan,
                'Active channel (1=Y, 0=N)': int(is_active),
                'Channel type (mesiotemporal or neocortical)': loc_group,
                'Channel type (SOZ or non-SOZ)': soz_group,
                'Pre-stimulus HFO no. on contact': hfo_pre_counts.get(chan, 0),
                'Arousal (/event) HFO no. on contact': hfo_post_counts.get(chan, 0)
            })

# --- 4. DATA SYNTHESIS AND NATURAL SORTING ---
df_hfo_final = pd.DataFrame(all_results)

# Apply natural sorting across Patient IDs, Event Types, and Channel Labels
df_hfo_final = df_hfo_final.sort_values(
    ['Patient ID', 'Event type', 'Event number', 'SEEG channel'],
    key=lambda col: col.str.replace('P', '').astype(int) if col.name == 'Patient ID' else col
)

