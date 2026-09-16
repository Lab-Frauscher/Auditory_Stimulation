# IED, HFO, Slow-Wave, and Spectral Power Ratio Analyses for Sleep Fragmentation Study

This repository provides the Python codebase and analysis pipeline used for our manuscript:
> **"Sleep fragmentation drives local, network-specific epileptic activity in human epilepsy"**

The pipeline processes stereo-EEG (SEEG) and scalp EEG recordings from patients with focal epilepsy to evaluate interictal epileptiform discharges (IEDs), high-frequency oscillations (HFOs), slow-wave co-occurrences, and spectral power ratios across auditory stimulation events and sleep arousals.

**Preprint:** [https://www.biorxiv.org/content/10.64898/2025.11.30.691386v1](https://www.biorxiv.org/content/10.64898/2025.11.30.691386v1)

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [System Requirements & Prerequisites](#system-requirements--prerequisites)
- [Installation Guide](#installation-guide)
- [Expected Runtimes](#expected-runtimes)
- [Reproducing the Analyses](#reproducing-the-analyses)
  - [Step 0: Automated Detection Pipeline (Optional)](#step-0-automated-detection-pipeline-optional)
  - [Step 1: Scalp EEG Power Ratios](#step-1-scalp-eeg-power-ratios)
  - [Step 2: IED and Slow-Wave Integration Analysis](#step-2-ied-and-slow-wave-integration-analysis)
  - [Step 3: High-Frequency Oscillations (HFO) Analysis](#step-3-high-frequency-oscillations-hfo-analysis)
  - [Step 4: IED Propagation & Thalamic Spectral Dynamics](#step-4-ied-propagation--thalamic-spectral-dynamics)
  - [Step 5: Generate Manuscript Figures](#step-5-generate-manuscript-figures)
- [External Detector Pipeline](#external-detector-pipeline)
- [Running on Custom Data](#running-on-custom-data)
- [References](#references)

---

## Repository Structure

```bash
├── LICENCE
├── README.md
├── active_channels.pkl
├── arousal_annotations.xlsx
├── channel_labels.xlsx
├── demo_data/
│   ├── auditory_stimulation_P11_demo.edf
│   └── auditory_stimulation_P14_demo.edf
├── detect_ied_hfo.py
├── detect_slow_wave.py
├── detector_results/
│   ├── p11_hfo.pkl
│   ├── p11_ied.pkl
│   ├── p11_slow-wave.pkl
│   ├── p14_hfo.pkl
│   ├── p14_ied.pkl
│   └── p14_slow-wave.pkl
├── figures/
│   ├── Figure_8a.pdf
│   ├── Figure_8b.pdf
│   ├── Figure_8c.pdf
│   ├── Figure_S16.pdf
│   ├── figures_8abc_S16.py
│   ├── Hannan_SourceData.xlsx
│   └── Hannan_Supplementary_SourceData.xlsx
├── hfo_analysis.py
├── ied_propagation_analysis.py
├── ied_slow-wave_power-ratio_analysis.py
├── powerratio_scalp_stimulation.pkl
├── requirements
├── spectral_power_scalp.py
├── thalamus_analysis.py
└── thalamus_channels.xlsx

```

---

## System Requirements & Prerequisites

### Operating System & Hardware

* **Operating System:** Linux, macOS, or Windows
* **RAM:** 8 GB minimum (16 GB recommended for parallel signal processing)
* **CPU:** Multi-core processor (recommended for multi-process detection in `detect_ied_hfo.py`)

### Python Dependencies

* **Python Version:** 3.10 or higher
* **Core Package Specifications:**
* `numpy` == 1.26.4
* `numba` ≥ 0.58.0
* `epycom` == 0.3.0
* `slow_wave_detector`
* `pyedflib`
* `scikit-learn`
* `scipy`
* `pandas`
* `openpyxl`
* `matplotlib`
* `seaborn`



---

## Installation Guide

> **Prerequisites:** 
> - Raw EEG demo recordings (`.edf`) use **Git Large File Storage (Git LFS)**. Ensure Git LFS is installed prior to cloning.

### 1. Clone the Main Repository
```bash
git clone [https://github.com/Lab-Frauscher/Auditory_Stimulation.git](https://github.com/Lab-Frauscher/Auditory_Stimulation.git)
cd Auditory_Stimulation
git lfs pull

```

### 2. Set Up Virtual Environment & Core Dependencies

```bash
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate

# Install core packages and fixed epycom/numba environment
pip install -r requirements.txt

```

### 3. Install External Slow-Wave Detector (GitLab SSH Required)

Clone and install the slow-wave detector package into your active environment:

```bash
# Clone slow-wave-detector repository
git clone git@gitlab.com:anphy_duke/lab-tools/slow-wave-detector.git

# Install slow-wave-detector in editable mode (or via its requirements)
cd slow-wave-detector
pip install -e .
cd ..

```

*Typical installation time: < 10 minutes.*

---

## Expected Runtimes

Evaluated on a standard desktop computer (Apple M-series / Intel i7 CPU, 16 GB RAM):

| Script | Approx. Runtime |
| --- | --- |
| `detect_ied_hfo.py` | ~346 s (~6 min) |
| `detect_slow_wave.py` | ~84 s (~1.5 min) |
| `spectral_power_scalp.py` | ~1 s |
| `ied_slow-wave_power-ratio_analysis.py` | ~1 s |
| `hfo_analysis.py` | ~2 s |
| `ied_propagation_analysis.py` | ~1 s |
| `thalamus_analysis.py` | ~1 s |
| `figures/figures_8abc_S16.py` | ~3 s |

**Total execution time:** ~7 minutes for full pipeline execution from raw EDFs (or < 10 seconds when using pre-computed detection outputs in `detector_results/`).

---

## Reproducing the Analyses

### Step 0: Automated Detection Pipeline (Optional)

*Pre-computed detector outputs for demo datasets are included in `detector_results/`. You can skip this step unless re-running detections directly from raw `.edf` recordings.*

To execute automated signal detection algorithms on demo files (`P11`, `P14`):

```bash
python detect_ied_hfo.py
python detect_slow_wave.py

```

---

### Step 1: Scalp EEG Power Ratios

Compute baseline pre-stimulus delta power ratio (0.5–4.0 Hz / 4.5–30.0 Hz) for scalp channels:

```bash
python spectral_power_scalp.py

```

---

### Step 2: IED and Slow-Wave Integration Analysis

Extract channel-level raw IED spike counts within 3-second baseline and post-event windows. Integrate event-level slow-wave presence, time-aligned scalp power ratios, sleep stages, SOZ designations, and anatomical classifications (mesiotemporal vs. neocortical):

```bash
python ied_slow-wave_power-ratio_analysis.py

```

---

### Step 3: High-Frequency Oscillations (HFO) Analysis

Extract ripple and fast-ripple events (≥ 80 Hz) from pre-computed detections, and quantify pre- and post-event HFO rates per channel:

```bash
python hfo_analysis.py

```

---

### Step 4: IED Propagation & Thalamic Spectral Dynamics

Group consecutive temporal spikes into discrete propagation events, and analyze thalamic SEEG power ratio dynamics:

```bash
python ied_propagation_analysis.py
python thalamus_analysis.py

```

---

### Step 5: Generate Manuscript Figures

Generate Figures 8a, 8b, 8c, and Supplementary Figure S16:

```bash
cd figures
python figures_8abc_S16.py

```

---

## External Detector Pipeline

The automated detection algorithms for IEDs, HFOs, and slow waves utilize established open-source signal processing frameworks:

* **IED Detection:** Detector methodology based on signal envelope distribution modeling ([Janca et al., 2015](https://doi.org/10.1007/s10548-014-0379-1?utm_source=gemini)). Implemented via `epycom.event_detection.spike.janca_detector`.
* **HFO Detection:** Automated detection of ripple and fast-ripple oscillations (≥ 80 Hz) ([von Ellenrieder et al., 2016](https://doi.org/10.1111/epi.13380?utm_source=gemini)). Implemented via `epycom.event_detection.hfo.nicolas_detector`.
* **Slow Wave Detection:** Automated detection of NREM slow oscillations ([Frauscher et al., 2015](https://doi.org/10.1093/brain/awv073?utm_source=gemini)). Implemented via `slow_wave_detector`.

---

## Running on Custom Data

To apply this pipeline to custom datasets:

1. **EDF Signal Preparation:** Place continuous patient recordings in `.edf` format under `demo_data/` following the naming convention `auditory_stimulation_P<ID>_demo.edf`.
2. **Event Annotations:** Update `arousal_annotations.xlsx` with UTC start timestamps (`utc_stim_start`, `utc_arousal_start`), event conditions (`evoked`, `spontaneous`, `control`), and sleep stages (`stage`, `stage_post`).
3. **Channel Mapping:** Define bipolar electrode contact pairs, clinical classifications (`SOZ` vs. `non-SOZ`), and anatomical locations (`mesiotemporal` vs. `neocortical`) in `channel_labels.xlsx` (and `thalamus_channels.xlsx` for subcortical channels).
4. **Detection Results:** Place external or recomputed detector outputs in `detector_results/` (`p<id>_ied.pkl`, `p<id>_hfo.pkl`, `p<id>_slow-wave.pkl`).

---


## References

1. Janca, R. et al. Detection of Interictal Epileptiform Discharges Using Signal Envelope Distribution Modelling: Application to Epileptic and Non-Epileptic Intracranial Recordings. *Brain Topogr* 28, 172–183 (2015). [https://doi.org/10.1007/s10548-014-0379-1](https://doi.org/10.1007/s10548-014-0379-1?utm_source=gemini)
2. von Ellenrieder, N., Frauscher, B., Dubeau, F. & Gotman, J. Interaction with slow waves during sleep improves discrimination of physiologic and pathologic high-frequency oscillations (80-500 Hz). *Epilepsia* 57, 869–878 (2016). [https://doi.org/10.1111/epi.13380](https://doi.org/10.1111/epi.13380?utm_source=gemini)
3. Frauscher, B. et al. Facilitation of epileptic activity during sleep is mediated by high amplitude slow waves. *Brain* 138, 1629–1641 (2015). [https://doi.org/10.1093/brain/awv073](https://doi.org/10.1093/brain/awv073?utm_source=gemini)
