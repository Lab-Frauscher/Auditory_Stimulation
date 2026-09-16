"""
ied_propagation_analysis.py

Description:
==================
This script implements a spatial feature extraction and event-clustering pipeline 
designed to evaluate interictal epileptiform discharge (IED) propagation dynamics
during sleep arousals. The analysis processes stereo-EEG 
(SEEG) recordings from patients with focal epilepsy, grouping detected temporal spikes 
into discrete propagation events and mapping their origins to specific anatomical regions 
(mesiotemporal vs. neocortical) and clinical zones (Seizure Onset Zone, SOZ vs. non-SOZ). 
Metrics are computed over paired 3-second baseline and post-arousal epochs across evoked, 
spontaneous, and control experimental events.

Key Processing Steps:
-------------------
1. Loading arousal event annotations, SEEG contact metadata, and detected IED peak timestamps.
2. Parsing bipolar channel structures and verifying anatomical/clinical concordance of constituent contacts.
3. Temporal clustering of individual spikes into unified IED propagation events using a 120 ms time-gap threshold.
4. Identifying event origin channels and measuring spatial propagation extent (channel participation count).
5. Window segmentation for baseline/pre-stimulus (-3 to 0 s) and post-arousal/post-stimulus (+0 to +3 s) conditions.
6. Calculating category-specific Global Spike Index (GSI) and mean spatial propagation metrics.
7. Partitioning aggregated outcomes by experimental condition (evoked, control, spontaneous) for downstream statistical analysis.

Dependencies:
-------------------
- numpy
- pandas
- openpyxl
"""

import pandas as pd
import numpy as np

# --- 1. DATA LOADING ---
# Load arousal event annotations and channel metadata
arousals = pd.read_excel('arousal_annotations.xlsx')
channel_labels = pd.read_excel('channel_labels.xlsx')

# --- 2. HELPER FUNCTIONS ---
def analyze_spikes(spikes_df, win_start, win_end, time_col='time_peak', chan_col='channel', min_gap_sec=0.12):
    """
    Groups individual spikes within a time window into discrete IED events based on a temporal 
    threshold and computes spatial propagation and initial channel origins.

    Parameters:
    -----------
    spikes_df : pd.DataFrame
        DataFrame containing detected spike timestamps and channel names.
    win_start : float
        Start time of the analysis window (in seconds).
    win_end : float
        End time of the analysis window (in seconds).
    time_col : str, default='time_peak'
        Column name for spike peak timestamps.
    chan_col : str, default='channel'
        Column name for recorded channel identifiers.
    min_gap_sec : float, default=0.12
        Maximum inter-spike interval threshold (120 ms) to group consecutive spikes 
        into the same propagation event.

    Returns:
    --------
    propagations : list of int
        Number of unique channels involved in each detected IED event.
    first_channels : list of str
        The initial channel identifier where each IED event originated.
    """
    # Filter spikes within the specified time window and sort chronologically
    mask = (spikes_df[time_col] >= win_start) & (spikes_df[time_col] <= win_end)
    win_spikes = spikes_df[mask].sort_values(by=time_col)
    
    if win_spikes.empty:
        return [], []
        
    propagations = []
    first_channels = []
    
    times = win_spikes[time_col].values
    chans = win_spikes[chan_col].values
    
    # Initialize the first IED event
    current_event_start = times[0]
    current_first_chan = chans[0]
    current_channels = {chans[0]}
    
    for i in range(1, len(times)):
        t = times[i]
        ch = chans[i]
        
        # If the gap exceeds the temporal threshold, finalize current event and start a new one
        if t > current_event_start + min_gap_sec:
            propagations.append(len(current_channels))
            first_channels.append(current_first_chan)
            
            current_event_start = t
            current_first_chan = ch
            current_channels = {ch}
        else:
            current_channels.add(ch)
            
    # Append the final event in the window
    propagations.append(len(current_channels))
    first_channels.append(current_first_chan)
    
    return propagations, first_channels

def get_bipolar_groups(bipolar_chan, loc_dict, soz_dict):
    """
    Parses a bipolar channel label and verifies whether both constituent unipolar 
    contacts share the same anatomical or clinical classification.

    Parameters:
    -----------
    bipolar_chan : str
        Bipolar channel string (e.g., 'A1-A2').
    loc_dict : dict
        Mapping of contact names to anatomical locations ('mesiotemporal' vs 'neocortical').
    soz_dict : dict
        Mapping of contact names to clinical zones ('SOZ' vs 'non-SOZ').

    Returns:
    --------
    groups : list of str
        List of classification groups ('mesiotemporal', 'neocortical', 'SOZ', 'non-SOZ') 
        that both contacts concordantly belong to.
    """
    groups = []

    # Extract individual unipolar contact labels from bipolar pair
    c1, c2 = [p.strip() for p in str(bipolar_chan).split('-')]
    
    # Retrieve metadata attributes (case-insensitive fallback)
    loc1, loc2 = str(loc_dict.get(c1, '')).lower(), str(loc_dict.get(c2, '')).lower()
    soz1, soz2 = str(soz_dict.get(c1, '')).lower(), str(soz_dict.get(c2, '')).lower()
    
    # Anatomical region concordance check
    if loc1 == 'mesiotemporal' and loc2 == 'mesiotemporal':
        groups.append('mesiotemporal')
    elif loc1 == 'neocortical' and loc2 == 'neocortical': 
        groups.append('neocortical')
        
    # Seizure Onset Zone (SOZ) concordance check
    if soz1 == 'soz' and soz2 == 'soz':
        groups.append('SOZ')
    elif soz1 == 'non-soz' and soz2 == 'non-soz':
        groups.append('non-SOZ')
        
    return groups

def calculate_group_metrics(first_channels, propagations, bipolar_group_map, group_totals, prefix=""):
    """
    Calculates Global Spike Index (GSI; spike rate normalized by total channel count) 
    and mean propagation extent for each anatomical/clinical category within a given window.

    Parameters:
    -----------
    first_channels : list of str
        List of originating channels for IED events in the window.
    propagations : list of int
        List of channel participation counts for IED events in the window.
    bipolar_group_map : dict
        Mapping of bipolar channels to their assigned group classifications.
    group_totals : dict
        Total count of valid bipolar channels present in each group for the given patient.
    prefix : str, default=""
        Prefix label added to output dictionary keys (e.g., 'Pre ' or 'Post ').

    Returns:
    --------
    metrics : dict
        Dictionary containing calculated GSI and average propagation per category.
    """
    counts = {'mesiotemporal': 0, 'neocortical': 0, 'SOZ': 0, 'non-SOZ': 0}
    props = {'mesiotemporal': [], 'neocortical': [], 'SOZ': [], 'non-SOZ': []}
    
    # Map detected events to functional/anatomical groups based on origin channel
    for ch, prop in zip(first_channels, propagations):
        grps = bipolar_group_map.get(ch, [])
        for g in grps:
            counts[g] += 1
            props[g].append(prop)
            
    # Calculate GSI (event count / channel count) and mean spatial propagation
    metrics = {}
    for g in counts.keys():
        total_ch = group_totals[g]
        gsi = counts[g] / total_ch if total_ch > 0 else np.nan
        mean_prop = np.mean(props[g]) if len(props[g]) > 0 else np.nan
        
        metrics[f'{prefix}{g} GSI'] = gsi
        metrics[f'{prefix}{g} Prop'] = mean_prop
        
    return metrics

# --- 3. DATA PROCESSING PIPELINE ---
all_results = []

for patient_id in ['P11', 'P14']:
    # Load automated IED detection results for the patient
    IED_id = pd.read_pickle('detector_results/' + patient_id.lower() + '_ied.pkl')
    
    # Filter and order patient-specific arousal events
    arousal_id = arousals.loc[arousals['ID'] == patient_id].sort_values(by=['utc_stim_start', 'utc_arousal_start'])
    
    # Restrict analysis to events within a 4-hour post-stimulation protocol window
    fragmentation_stop = arousal_id['utc_stim_start'].max()
    arousal_id = arousal_id.loc[(arousal_id['utc_arousal_start'] <= fragmentation_stop + 4 * 3600) | (arousal_id['utc_arousal_start'].isna())]
    
    # Filter channel labels for current patient
    channel_label_id = channel_labels.loc[channel_labels['Patient ID'] == patient_id]

    # --- PATIENT-SPECIFIC CHANNEL MAP PREPARATION ---
    # Construct contact lookup dictionaries
    loc_dict = dict(zip(channel_label_id['Contact'].astype(str).str.upper().str.strip(), 
                        channel_label_id['mesiotemporal vs neocortical'].astype(str).str.strip()))
    soz_dict = dict(zip(channel_label_id['Contact'].astype(str).str.upper().str.strip(), 
                        channel_label_id['SOZ vs non-SOZ'].astype(str).str.strip()))
    
    # Identify unique recorded bipolar channels for the current dataset
    unique_bipolar_chans = IED_id['channel'].unique()
    
    # Map bipolar channels to group categories and compute channel totals per category
    bipolar_group_map = {}
    group_channel_totals = {'mesiotemporal': 0, 'neocortical': 0, 'SOZ': 0, 'non-SOZ': 0}
    
    for chan in unique_bipolar_chans:
        grps = get_bipolar_groups(chan, loc_dict, soz_dict)
        bipolar_group_map[chan] = grps
        for g in grps:
            group_channel_totals[g] += 1
            
    total_channels = len(unique_bipolar_chans) 
    
    patient_spikes_df = IED_id[['time_peak', 'channel']]
    
    # --- AROUSAL EVENT PROCESSING ---
    for index, row in arousal_id.iterrows():
        a_type = row['type']
        
        stim_time = row['utc_stim_start']
        arousal_start = row['utc_arousal_start']
        
        # Window segmentation logic based on arousal condition
        if a_type == 'evoked':
            # Pre-stimulus (3s baseline) vs. Post-arousal onset (3s window)
            prop_pre, first_ch_pre = analyze_spikes(patient_spikes_df, stim_time - 3, stim_time)
            prop_post, first_ch_post = analyze_spikes(patient_spikes_df, arousal_start, arousal_start + 3)
            
        elif a_type == 'control':
            # Pre-stimulus baseline vs. Time-matched post-stimulus control window (offset by 1.076s)
            prop_pre, first_ch_pre = analyze_spikes(patient_spikes_df, stim_time - 3, stim_time)
            win_start_post = stim_time + 1.076
            prop_post, first_ch_post = analyze_spikes(patient_spikes_df, win_start_post, win_start_post + 3)
            
        elif a_type == 'spontaneous':
            # Pre-arousal (3s baseline) vs. Post-arousal onset (3s window)
            prop_pre, first_ch_pre = analyze_spikes(patient_spikes_df, arousal_start - 3, arousal_start)
            prop_post, first_ch_post = analyze_spikes(patient_spikes_df, arousal_start, arousal_start + 3)
            
        else:
            continue
            
        # Compute category-level metrics (GSI and spatial propagation)
        pre_group_metrics = calculate_group_metrics(first_ch_pre, prop_pre, bipolar_group_map, group_channel_totals, prefix="Pre ")
        post_group_metrics = calculate_group_metrics(first_ch_post, prop_post, bipolar_group_map, group_channel_totals, prefix="Post ")
        
        # Compile result row for summary dataset
        result_row = {
            'Patient ID': patient_id,
            'Type': a_type,
            'Stimulus name': row.get('stim_intensity', np.nan),
            'Pre-stimulus Sleep Stage': row.get('stage', np.nan),
            'Post-stimulus Sleep Stage': row.get('stage_post', np.nan),
            
            # Global overall metrics
            'Pre-stimulus IED Propagation Channels': np.mean(prop_pre) if len(prop_pre) > 0 else np.nan,
            'Arousal IED Propagation Channels': np.mean(prop_post) if len(prop_post) > 0 else np.nan,
            
            # Subgroup spatial propagation metrics
            'Pre-stimulus mesiotemporal IED Propagation Channels': pre_group_metrics['Pre mesiotemporal Prop'],
            'Arousal mesiotemporal IED Propagation Channels': post_group_metrics['Post mesiotemporal Prop'],
            
            'Pre-stimulus neocortical IED Propagation Channels': pre_group_metrics['Pre neocortical Prop'],
            'Arousal neocortical IED Propagation Channels': post_group_metrics['Post neocortical Prop'],
            
            'Pre-stimulus SOZ IED Propagation Channels': pre_group_metrics['Pre SOZ Prop'],
            'Arousal SOZ IED Propagation Channels': post_group_metrics['Post SOZ Prop'],
            
            'Pre-stimulus non-SOZ IED Propagation Channels': pre_group_metrics['Pre non-SOZ Prop'],
            'Arousal non-SOZ IED Propagation Channels': post_group_metrics['Post non-SOZ Prop'],
        }
        
        all_results.append(result_row)

# --- 4. DATA SYNTHESIS AND SEPARATION ---
df_results = pd.DataFrame(all_results)

# Partition dataset into distinct experimental conditions for statistical analysis
evoked_df = df_results[df_results['Type'] == 'evoked'].copy()
control_df = df_results[df_results['Type'] == 'control'].copy()
spontaneous_df = df_results[df_results['Type'] == 'spontaneous'].copy()
