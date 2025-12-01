# IED analyses for sleep fragmentation study

This repository provides the code used to reproduce the main analyses for our manuscript *‘Sleep fragmentation drives local, network-specific epileptic activity in human epilepsy’*. This code included here is used for analysing interictal epileptiform discharges (IEDs) and thalamic power ratio from SEEG data in relation to auditory stimulation.

You can find our preprint here:

<http://insert_url>

## Table of Contents

The directory structure is as follows:

-   Folders
-   Prerequisites
-   Installation
-   Expected runtime
-   Reproducing the analyses
-   Running the code on your own data
-   References

## Folders

``` bash
├── LICENSE
├── README.md
├── arousal_annotations.xlsx
├── channel_labels.xlsx
├── demo_data
│   ├── P11.mat
│   └── P14.mat
├── figures
│   ├── Figure_7a.pdf
│   ├── Figure_7b.pdf
│   ├── Figure_7c.pdf
│   ├── Figure_S2.pdf
│   ├── figures_7abc_S2.py
│   ├── Hannan_SourceData.xlsx
│   └── Hannan_Supplementary_SourceData.xlsx
├── functions 
│   ├── postprocessing_v3.m
│   └── spike_detector_hilbert_v25.m
├── IED_analysis.m
├── thalamus_channels.xlsx
└── thalamus_powerratio.py
```

## Prerequisites

Operating system: Windows, macOS, or Linux

### MATLAB

-   MATLAB R2022b or higher
-   Required toolboxes:

1.  Signal Processing Toolbox
2.  Statistics and Machine Learning Toolbox

### Python

-   Python 3.10.11 or higher
-   Required libraries:

1.  numpy ≥ 1.25.0
2.  pandas ≥ 2.2.3
3.  scipy ≥ 1.11.1
4.  h5py ≥ 3.9.0
5.  matplotlib ≥ 3.7.2
6.  seaborn ≥ 0.13.2

## Installation

Download the repository using the following *git* command:

*git clone <https://github.com/Lab-Frauscher/Auditory_Stimulation>*

### MATLAB setup:

Install each prerequisite toolbox by first opening MATLAB R2022b and navigating to Home \> Add-ons \> Get Add-ons. Type the name of each respective toolbox into the search bar and install it.

### Python setup:

It is strongly recommended to create and activate a virtual environment first. Then, install all dependencies for Python using *pip*:

*pip install -r requirements.txt*

Typical installation time: 2 hours

## Expected runtime

Time to run *IED_analysis.m*: 101.6983 seconds

Time to run *thalamus_powerratio.py*: 266.6811 seconds

Time to run *figures_7abc_S2.py*: 3.5569 seconds

## Reproducing the analyses

The demo scripts use sample SEEG data from two patients (*P11.mat*, *P14.mat*) (\~20 min sampled at 2048 Hz), which can be found in the demo_data subdirectory.

The structure of the *Pxx.mat* file is as follows:

-   *bp_data*: the core bipolar SEEG signal data, stored as a double array. The dimensions are (number of channels) x (number of samples).
-   *bp_ch_names*: a cell array providing the labels for each bipolar channel.
-   *fs*: a double containing the signal's sampling frequency in Hz.
-   *n_bp*: a double specifying the total number of bipolar channels.
-   *n_samples2*: a double specifying the total number of recorded samples.

Set MATLAB workspace to the repository folder, and get started by running the demo script for IED analysis on the sample data:

*IED_analysis.m*

The script creates variables *arousals_GIED* and *controls_GIED* (normalized spike rates for pre-stimulus and arousal/post-stimulus windows), and *arousals_propagation* and *controls_propagation* (propagation results).

Next, open a Python environment with installed requirements and run the demo script for thalamus power ratio analysis

*thalamus_powerratio.py*

This script uses the same sample data but focuses only on stimuli that induced arousals. The variable *power_ratio_arousals* will show a relative change in the power ratio.

If you want to reproduce the result figures in the manuscript, you need to go to the folder

``` bash
figures
├── Figure_7a.pdf
├── Figure_7b.pdf
├── Figure_7c.pdf
├── Figure_S2.pdf
├── figures_7abc_S2.py
├── Hannan_SourceData.xlsx
└── Hannan_Supplementary_SourceData.xlsx
```

If you want to reproduce the result figures in the manuscript, the necessary results applied to all patients are in the *Hannan_SourceData.xlsx* and *Hannan_Supplementary_SourceData.xlsx* files in the figures folder. Run the Python script to generate selected figures:

*figures_7abc_S2.py*

Note:

The function *spike_detector_hilbert_v25.m* for IED detection was previously developed for earlier work (<https://github.com/Lab-Frauscher/Spike-Gamma>) and is based on the method described by *Janca et al., 2015*. The function *postprocessing_v3.m* is an updated version of earlier code, adapted for the analyses reported here. All other scripts were written specifically for this study.

## Running the code on your own data

To use these scripts with your own SEEG data:

-   Extract Data: extract the specific segment of interest from your SEEG data that corresponds to the auditory stimulation protocol.
-   Apply Montage: apply the bipolar montage to the extracted signal.
-   Annotate Arousals: annotate or label all arousals present within your processed data segment.
-   Determine Anatomy: determine the anatomical position for each channel, including identifying channels in the thalamus and distinguishing between mesiotemporal and neocortical channels.
-   Separate Channels: separate the channels into groups Seizure Onset Zone (SOZ) and non-SOZ channels, based on clinical findings.
-   Run scripts *IED_analysis.m* and *thalamus_powerratio.py*.

## References

Janca, R., Jezdik, P., Cmejla, R. et al. Detection of Interictal Epileptiform Discharges Using Signal Envelope Distribution Modelling: Application to Epileptic and Non-Epileptic Intracranial Recordings. Brain Topogr 28, 172--183 (2015). <https://doi.org/10.1007/s10548-014-0379-1>
