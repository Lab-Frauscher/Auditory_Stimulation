"""
figures_8abc_S16.py

Description:
==================
This script plots figures to visualize the modulation of Interictal Epileptiform
Discharges (IEDs) in human iEEG data. It loads pre-aggregated and normalized IED
count data from Excel source files and generates four plots showing IED counts
binned to 500ms intervals in time.

Key Processing Steps:
-------------------
1. Loading the necessary IED count data for cohort and individual
patient analysis.
2. Binning to confirm the use of a 500 ms bin length for all time-relative measurements.
3. Figures 8a-8c (Group Plots): The normalized IED count for the main cohort across three conditions:
   Figure 8a: IED count aligned to the arousal onset.
   Figure 8b: IED count aligned to the stimulation onset (for stimulations that
              resulted in arousal), overlaid with a line showing the percentage
              of active arousals.
   Figure 8c: IED count aligned to the stimulation onset for control stimulations.
4. Supplementary Figure 16: IED count aligned to arousal onset for each individual patient ID

Dependencies:
-------------------
- numpy
- pandas
- matplotlib.pyplot
- seaborn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Data Reading
file = 'Hannan_SourceData.xlsx'

fig_7a = pd.read_excel(file, sheet_name='Figure 8', header=4, usecols='B:D').rename(columns={
    'time relative to arousal onset (s, 500-ms bin centre)': 'time relative to arousal onset (s)'}).dropna(how='all')

fig_7b = pd.read_excel(file, sheet_name='Figure 8', header=4, usecols='H:J').rename(columns={
    'time relative to stimulus onset (s, 500-ms bin centre)': 'time relative to stimulus onset (s)',
    'IED count (normalised).1': 'IED count (normalised)'}).dropna(how='all')
fig_7b_line = pd.read_excel(file, sheet_name='Figure 8', header=4, usecols='F:G').dropna(how='all')

fig_7c = pd.read_excel(file, sheet_name='Figure 8', header=4, usecols='L:N').rename(columns={
    'time relative to stimulus onset (s, 500-ms bin centre).1': 'time relative to stimulus onset (s)',
    'IED count (normalised).2': 'IED count (normalised)'}).dropna(how='all')

file_supp = 'Hannan_Supplementary_SourceData.xlsx'
fig_s2 = pd.read_excel(file_supp, sheet_name='Supplementary Figure 16', header=3, usecols='B:D').rename(columns={
    'time relative to arousal onset (s, 500-ms bin centre)': 'time relative to arousal onset (s)'})

# bin length 500 ms
bin_len = 0.5

# Figure 8a: Bar Plot of IED Count During Arousals Relative to Arousal Onset
g = sns.catplot(
    data=fig_7a,
    kind="bar",
    x="time relative to arousal onset (s)",
    y="IED count (normalised)",
    color="#1151b8",
    height=2.7,
    aspect=1.8,
    estimator=np.mean,
    errorbar='se',
    capsize=0.15,
    errwidth=1.5,
    errcolor='black',
    sharex=True,
    sharey=False
)

g.tick_params(axis='both', length=0)

ax = g.axes.flatten()[0]
for ax in g.axes.flatten():
    ax.set_xticks(np.arange(0, 11/bin_len, 1/bin_len)-bin_len*0.8,np.arange(-5, 6, 1))
    ax.set_xlim(-1, 10/bin_len)
    ax.set_ylim(0.4, 1.1)
    
    for patch in ax.patches:
        val = patch.get_height()
        patch.set_y(0.4)
        patch.set_height(val - 0.4)
    

plt.tight_layout()
plt.show()

# Figure 8b: Bar Plot of IED Count During Arousals Relative to Stimulus Onset
g = sns.catplot(
    data=fig_7b,
    kind="bar",
    x="time relative to stimulus onset (s)",
    y="IED count (normalised)",
    color="#1151b8",
    height=2.7,
    aspect=1.8,
    estimator=np.mean,
    errorbar='se',
    capsize=0.15,
    errwidth=1.5,
    errcolor='black',
    sharex=True,
    sharey=False
)

g.tick_params(axis='both', length=0)

ax = g.axes.flatten()[0]

for ax in g.axes.flatten():
    ax.set_xticks(np.arange(0, 12/bin_len, 1/bin_len)-bin_len*0.8,np.arange(-3, 9, 1))
    ax.set_xlim(-1, 11/bin_len)
    ax.set_ylim(0.4, 1.1)
    
    for patch in ax.patches:
        val = patch.get_height()
        patch.set_y(0.4)
        patch.set_height(val - 0.4)
    
ax2 = ax.twinx()
ax2.plot(((fig_7b_line['time relative to  stimulus onset (s)']+3)/bin_len)-bin_len*0.9, fig_7b_line['% active arousals'], color="black", linewidth=1.5)
ax2.set_ylabel("% active arousals")
ax2.tick_params(axis='both', length=0)
ax2.set_ylim(0, 140)
ax2.spines['top'].set_visible(False)

plt.tight_layout()
plt.show()

# Figure 8c: Bar Plot of IED Count for Controls Relative to Stimulus Onset
g = sns.catplot(
    data=fig_7c,
    kind="bar",
    x="time relative to stimulus onset (s)",
    y="IED count (normalised)",
    color="#1151b8",
    height=2.7,
    aspect=1.8,
    estimator=np.mean,
    errorbar='se',
    capsize=0.15,
    errwidth=1.5,
    errcolor='black',
    sharex=True,
    sharey=False
)

ax = g.axes.flatten()[0]

g.tick_params(axis='both', length=0)

for ax in g.axes.flatten():
    ax.set_xticks(np.arange(0, 12/bin_len, 1/bin_len)-bin_len*0.8,np.arange(-3, 9, 1))
    ax.set_xlim(-1, 11/bin_len)
    ax.set_ylim(0.4, 1.1)
    
    for patch in ax.patches:
        val = patch.get_height()
        patch.set_y(0.4)
        patch.set_height(val - 0.4)
    
plt.tight_layout()
plt.show()

# Supplementary Figure 16: Bar Plot of IED Count During Arousals Relative to Arousal Onset for Separated Patiens
g = sns.catplot(
    data=fig_s2,
    kind="bar",
    x="time relative to arousal onset (s)",
    y="IED count (normalised)",
    row="patient ID",
    row_order=['P17', 'P16', 'P15', 'P14', 'P13', 'P12', 'P11', 'P10', 'P9', 'P8',
               'P7', 'P6', 'P5', 'P4','P3', 'P2', 'P1'][::-1],
    height=0.4,
    aspect=9,
    color="#1151b8",
    sharex=True,
    sharey=False
)

g.set_titles('')

g.tick_params(axis='both', length=0)

for ax, patient_id in zip(g.axes.flatten(), g.row_names):
    ax.set_ylabel(patient_id, rotation=0, labelpad=10, fontsize=10, va='center', ha='right')
    ax.set_yticks([])
    ax.set_xticks(np.arange(0, 11/bin_len, 1/bin_len)-bin_len*0.8,np.arange(-5, 6, 1))
    ax.set_xlim(-1, 10/bin_len)

g.fig.text(x=0.01, y=0.5, s="IED count (normalised)", va='center', rotation=90, fontsize=11)

plt.subplots_adjust(top=0.97, bottom=0.08, left=0.18)
plt.show()
