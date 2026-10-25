# Neurocomputational Analysis of Decision-Making Under Uncertainty

This repository contains the MATLAB analysis scripts used for the behavioral, computational, EEG time-frequency, mixed-effects, sensitivity, and reliability analyses reported in the project.

The analysis combines:

- Cumulative Prospect Theory (CPT) modeling
- trial-wise subjective value estimation
- behavioral mixed-effects modeling
- EEG time-frequency analysis in the beta band
- frontoparietal ROI analyses
- ambiguity-sensitivity analyses
- split-half reliability analysis
- publication figure generation

The main analyses focus on how model-derived decision conflict relates to response timing, how this relationship differs across cognitive-control groups, and whether frontoparietal beta modulation is associated with these differences.

## Repository structure

```text
.
├── TF_Analysis.m
├── BETA_FIGURES.m
├── Cumulative_Prospect_Theory.m
├── LMM1.m
├── LMM2.m
├── AMBIGUITY_PRIOR_SENSITIVITYn.m
├── HA_RP_RELIABILITY.m
└── README.md
```

## Analysis overview

### 1. EEG time-frequency and ROI extraction

**Script:** `MASTER_TF_ROI_GitHub_clean.m`

This script performs the main EEG time-frequency analysis from preprocessed EEGLAB `.set` files.

Main steps:

- loads the final cleaned EEG datasets
- identifies task-specific event markers
- epochs trials for:
  - Risky Positive (RP)
  - Risky Negative (RN)
  - Low Ambiguity (LA)
  - High Ambiguity (HA)
- applies a time-domain baseline correction from -200 to 0 ms
- runs EEGLAB `newtimef`
- extracts beta-band activity from 13-30 Hz
- averages beta power over the 0-1000 ms analysis window
- stores channel-wise beta values for each participant and condition
- calculates the predefined ROI9 and ROI6 summaries
- computes the participant-level HA-RP beta modulation index
- performs the paired ROI9 RP-versus-HA comparison

The primary ROI9 contains:

```text
F3, F4, Fz, FCz, P3, P4, Pz, CPz, POz
```

The ROI6 sensitivity set contains:

```text
F3, F4, Fz, P3, P4, Pz
```

Missing channels, failed time-frequency estimates, or insufficient trial counts remain `NaN`; they are not replaced with zero.

Main outputs include:

```text
TF_OUTPUT_FINAL.mat
ROI_4Cond_pt_beta_FINAL.csv
ROI9_ROI6_summary_FINAL.csv
trial_counts.csv
task_status.csv
problem_log.csv
```

---

### 2. Beta time-course and scalp figures

**Script:** `BETA_FIGURES_GitHub_clean.m`

This script generates the EEG figures from the saved time-frequency output.

It includes:

- RP and HA frontoparietal beta time courses
- HA-RP beta modulation time courses
- group-specific Low-CC and High-CC HA-RP waveforms
- absolute scalp beta topographies
- HA-RP scalp difference maps
- spatiotemporal scalp maps across:
  - 200-400 ms
  - 400-600 ms
  - 600-800 ms
  - 800-1000 ms
- common color scaling across topographies
- source-data export for the plotted results

A light smoothing/interpolation step is used only for waveform display. Statistical analyses and exported source values remain based on the original data.

---

### 3. CPT model figures

**Script:** `CPT_FIGURES_GitHub_clean.m`

This script generates figures describing the computational model and its empirical behavior.

It produces:

- a schematic of the CPT analysis pipeline
- the probability-weighting function
- the outcome-utility function
- empirical subjective-value distributions
- participant-level model performance
- risky and ambiguous prediction accuracy
- McFadden pseudo-R² for risky-trial model fit
- choice-probability calibration curves

The CPT model uses participant-specific parameters:

```text
alpha   outcome sensitivity
lambda  loss aversion
gamma   probability weighting
tau     choice stochasticity / temperature
```

For trial-level subjective value:

```text
SV = w(p_gain) * gain^alpha
     - lambda * w(p_loss) * loss^alpha
```

The absolute subjective value, `|SV|`, is used as a model-derived index of proximity to subjective indifference. Smaller `|SV|` corresponds to greater modeled decision conflict.

---

### 4. Behavioral mixed-effects analysis: LMM1

**Script:** `LMM1_GitHub_clean.m`

LMM1 tests whether the relationship between model-derived conflict and response timing differs between cognitive-control groups.

The fitted model is:

```text
Log_RT ~ Z_Abs_SV * CC_Group + (Z_Abs_SV | ID)
```

The model is fitted separately for risky and ambiguous trials.

The final extreme-groups sample contains:

```text
N = 38
19 High-CC participants
19 Low-CC participants
```

The main behavioral predictor is standardized absolute subjective value:

```text
Z_Abs_SV = z(|SV|)
```

For plotting only, the displayed conflict axis is reversed:

```text
Displayed_Conflict = -Z_Abs_SV
```

Thus, moving rightward on the displayed axis corresponds to greater model-derived conflict.

The script also produces:

- model-based prediction curves
- 95% confidence intervals
- participant-balanced observed bins
- condition-specific simple slopes
- High-versus-Low CC interaction contrasts
- publication-ready figures and source-data files

---

### 5. Neurobehavioral mixed-effects analysis: LMM2

**Script:** `LMM2_Q1_GitHub_clean.m`

LMM2 tests whether the conflict-response-time relationship varies as a function of frontoparietal beta modulation and cognitive-control group.

The primary model is:

```text
Log_RT ~ Z_Abs_SV * Z_Neural_Braking * CC_Group + (1 | ID)
```

where:

```text
Z_Neural_Braking
```

is the standardized participant-level HA-RP beta modulation index derived from the primary ROI9.

The analysis uses ambiguous trials from the same final N=38 extreme-groups sample.

The script includes:

- the primary three-way interaction
- gatekept within-group follow-up models
- ROI6 sensitivity analysis
- prediction-based fixed-effect and full-model R²
- model-based interaction plots
- within-group beta-by-conflict coefficients
- ROI9-versus-ROI6 robustness comparison
- source-data export

The figures use:

```text
Displayed_Conflict = -Z_Abs_SV
```

for interpretation, while the statistical model itself is fitted using `+Z_Abs_SV`.

---

### 6. Ambiguity-prior sensitivity analysis

**Script:** `AMBIGUITY_PRIOR_SENSITIVITY_GitHub_clean.m`

This script evaluates whether the main LMM1 and LMM2 results depend on assumptions used to compute subjective value on ambiguous trials.

CPT parameters are estimated once from risky trials only.

Ambiguous-trial subjective values are then recomputed under nine combinations of:

```text
Assumed gain probability:
0.40
0.50
0.60

Assumed hidden magnitude:
10
12.5
15
```

with:

```text
p_loss = 1 - p_gain
```

For each scenario, the script reruns:

```text
LMM1:
Log_RT ~ Z_Abs_SV * CC_Group + (Z_Abs_SV | ID)

LMM2:
Log_RT ~ Z_Abs_SV * Z_HA_RP_Beta * CC_Group + (1 | ID)
```

The primary specification is:

```text
p_gain = 0.50
p_loss = 0.50
hidden magnitude = 12.5
```

The sensitivity analysis records:

- regression coefficient
- standard error
- 95% confidence interval
- F statistic
- degrees of freedom
- p value
- number of trials
- whether the effect direction matches the primary analysis

Main output:

```text
AMBIGUITY_PRIOR_SENSITIVITY_RESULTS.csv
```

---

### 7. HA-RP split-half reliability

**Script:** `HA_RP_RELIABILITY_GitHub_clean.m`

This script estimates the reliability of the participant-level HA-RP beta modulation index.

The analysis reruns the same `newtimef` procedure used in the primary EEG pipeline on randomly split trial halves.

Reliability is evaluated for:

```text
HA-RP difference
HA alone
RP alone
```

and for:

```text
full EEG-valid sample: N = 57
extreme-groups sample: N = 38
```

The script first reconstructs the full-condition ROI9 values and checks them against the primary ROI output. Reliability is calculated only if the reconstructed values match the primary pipeline within the specified numerical tolerance.

For each repeated split:

1. RP trials are randomly divided into two balanced halves.
2. HA trials are independently divided into two balanced halves.
3. ROI9 beta is recomputed separately for each half.
4. half-to-half correlations are calculated across participants.
5. correlations are corrected using the Spearman-Brown formula.

The default number of random splits is:

```text
N_SPLITS = 50
```

Outputs include:

```text
HA_RP_reliability_summary_EXACT.csv
HA_RP_reliability_all_splits_EXACT.csv
HA_RP_full_condition_reconstruction_EXACT.csv
HA_RP_reliability_EXACT.mat
```

## Data processing conventions

Several preprocessing and quality-control decisions are shared across the scripts.

### Participant pool

The final EEG-valid pool contains 57 participants.

The primary extreme-groups analyses use 38 participants:

```text
19 High CC
19 Low CC
```

The middle cognitive-control tertile is excluded from the primary LMM1 and LMM2 extreme-group analyses.

### Reaction-time corrections

The final behavioral pipeline includes:

- removal of participant ID 9
- correction of confirmed decimal-shift RT errors below 1 s by multiplying them by 10
- removal of one confirmed erroneous ambiguous-trial RT for participant 67
- log transformation of corrected RT values

### Cognitive-control grouping

Cognitive-control groups are defined using the normalized Stroop interference score.

Larger Stroop interference values indicate greater interference and lower cognitive control.

The final extreme-group labels are kept identical across the main LMM1 and LMM2 analyses.

### Subjective value and conflict

The computational conflict predictor is based on the absolute CPT-derived subjective value:

```text
|SV|
```

Smaller values indicate greater proximity to subjective indifference.

The statistical models are fitted using standardized positive `|SV|`:

```text
Z_Abs_SV = z(|SV|)
```

For visualization only:

```text
Displayed_Conflict = -Z_Abs_SV
```

so that larger displayed values correspond to greater modeled conflict.

## Software requirements

The analysis was developed in MATLAB.

Required components include:

- MATLAB
- Statistics and Machine Learning Toolbox
- Optimization Toolbox
- EEGLAB
- EEGLAB `newtimef`
- standard EEGLAB dataset and epoching functions

Main EEGLAB functions used include:

```text
pop_loadset
eeg_checkset
pop_epoch
pop_rmbase
pop_select
newtimef
topoplot
```

## Running the analyses

The scripts are intended to be run approximately in the following order:

```text
1. MASTER_TF_ROI_GitHub_clean.m
2. BETA_FIGURES_GitHub_clean.m
3. CPT_FIGURES_GitHub_clean.m
4. LMM1_GitHub_clean.m
5. LMM2_Q1_GitHub_clean.m
6. AMBIGUITY_PRIOR_SENSITIVITY_GitHub_clean.m
7. HA_RP_RELIABILITY_GitHub_clean.m
```

Before running a script, update the file and output paths near the top of the file as needed.

Some scripts depend on intermediate files created by earlier stages.

## Reproducibility

The scripts include internal checks for:

- expected participant counts
- participant IDs
- trial counts
- ROI channel availability
- EEG channel matching
- time-frequency grid consistency
- missing values
- reconstruction of primary ROI metrics
- consistency of extreme-group membership

The split-half reliability analysis uses a fixed random seed to make repeated random splits reproducible.

## Notes

The code in this repository reflects the final analysis pipeline used for the project. Local file paths in the scripts are specific to the original analysis environment and should be replaced before running the code on another system.

Raw EEG and behavioral data are not included in this repository unless separately provided.

## Citation

If you use this code, please cite the associated manuscript once the final publication information is available.

Manuscript title:

> **Frontoparietal beta modulation is associated with cognitive-control differences in conflict-sensitive response timing under ambiguity**
