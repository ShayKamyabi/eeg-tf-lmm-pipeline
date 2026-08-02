# eeg-tf-lmm-pipeline
An analytical pipeline for stimulus-locked EEG data, including preprocessing, time-frequency analysis, and Linear Mixed Models (LMM)
# Neural Dynamics During Decision-Making: EEG Analysis Pipeline

This repository contains the MATLAB codebase and analytical framework for processing electrophysiological data (EEG) associated with our decision-making study. The pipeline encompasses everything from raw data preprocessing to final statistical modeling, designed to evaluate brain wave oscillations across experimental groups.

## Project Overview

The experimental methodology involves a decision-making task featuring both **single-option** and **two-option** conditions. The scripts provided here handle the full analysis pipeline for these conditions, including:
* Preprocessing and artifact removal using Independent Component Analysis (ICA).
* Time-frequency analysis of neural oscillations.
* Statistical evaluation using Linear Mixed Models (LMM) and mixed ANOVAs.

## Repository Structure

The code is organized sequentially to match the order of analysis:

* **`01_time_frequency/`**: Scripts to extract and analyze brain wave oscillations and spectral power across the chosen single-option trials.
* **`02_ult_LMM1/`**: 
* **`03_ult_LMM2/`**: The statistical modeling framework, including the setup and execution of the mixed ANOVAs and both Linear Mixed Models (`02-ult_lmm1` and `03-ult_lmm2`).

## Prerequisites and Requirements

To run this code, you will need the following software and toolboxes:
* **MATLAB** (Version R20XXx or higher recommended)
* **EEGLAB** (for preprocessing and ICA)
* MATLAB Statistics and Machine Learning Toolbox (for `fitlme` / LMM and ANOVA functions)

## Data Availability

**Note on Data:** Due to GitHub's file size limitations, the raw and preprocessed EEG datasets (`.set`, `.fdt`, or raw formats) are not hosted in this repository. 
* To run these scripts, you must place your data files in a local directory named `data/` at the root of this project, or adjust the file paths in the scripts accordingly. 
* [If applicable, insert link to where your data is hosted, e.g., OpenNeuro, OSF, or state "Data available upon reasonable request".]

## Usage 

1. Clone this repository to your local machine.
2. Ensure you have the required dependencies installed and added to your MATLAB path.
3. Place your raw datasets into your local working directory.
4. Run the scripts in numerical order, starting with `01_preprocessing/main_preprocess.m`. 

## License

[Standard academic text: e.g., This project is licensed under the MIT License - see the LICENSE file for details.]
