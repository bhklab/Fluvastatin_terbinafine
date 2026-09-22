# PDO Synergy Analysis

## Overview
This repository contains R scripts for analyzing drug combination synergy in patient-derived organoid (PDO) models using the SynergyFinder package. The analysis focuses on Fluvastatin (FLV) and Terbinafine (TBF) combinations across multiple breast cancer cell lines.

## Scripts

### SynergyFinder.R
Main script for drug combination synergy analysis using the SynergyFinder R package.

**Key Features:**
- Drug combination synergy calculation (ZIP, Bliss methods)
- Dose-response curve generation
- Sensitivity analysis across cell lines
- Synergy heatmap visualization
- Export of results and synergy scores

## Installation Requirements

### R Packages
```r
# Core packages
library(readr)
library(synergyfinder)
library(reticulate)

# Install SynergyFinder from Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("synergyfinder")
```

## Data Format

### Input Data Structure
The script expects a CSV file with the following columns:
- `date`: Experiment date
- `conc1`: Drug 1 concentration (FLV)
- `conc2`: Drug 2 concentration (TBF)
- `response`: Cell viability measurement
- `repeats`: Experimental replicate number
- `drug1`: Primary drug identifier
- `drug2`: Secondary drug identifier
- `conc_unit1`: Concentration unit for drug 1
- `conc_unit2`: Concentration unit for drug 2
- `block_id`: Cell line identifier

### Data Requirements
- **Drug combinations**: FLV + TBF pairs
- **Single drug controls**: FLV only (conc2 = 0)
- **Cell lines**: Multiple PDO models (excluding block_id 64)
- **Replicates**: Multiple experimental repeats per condition

## Usage

### Basic Workflow
1. **Load data**: Read FLV_TBF_data.csv
2. **Data preparation**: Combine combination and single drug data
3. **Reshape data**: Convert to SynergyFinder format
4. **Calculate synergy**: Apply ZIP and Bliss methods
5. **Generate plots**: Dose-response curves and synergy heatmaps
6. **Export results**: Save reports and synergy scores

### Key Functions
- `ReshapeData()`: Convert data to SynergyFinder format
- `CalculateSynergy()`: Compute synergy scores using multiple methods
- `CalculateSensitivity()`: Analyze drug sensitivity across cell lines
- `PlotDoseResponse()`: Generate dose-response visualizations
- `PlotSynergy()`: Create synergy heatmaps

## Output Files

### Generated Plots
- **Dose-response curves**: Individual and combined drug effects
- **Synergy heatmaps**: ZIP synergy scores across concentration ranges
- **Sensitivity plots**: Drug response patterns by cell line

### Data Exports
- **Synergy report**: Summary of synergy scores and parameters
- **Raw data**: Processed synergy scores for further analysis
- **Cell line comparisons**: Response patterns across different PDO models

## Analysis Methods

### Synergy Calculation
- **ZIP (Zero Interaction Potency)**: Evaluates drug combination effects
- **Bliss Independence**: Assesses additive vs. synergistic interactions
- **Baseline correction**: Non-corrected baseline approach
- **Iterations**: 10 iterations for robust calculations

### Sensitivity Analysis
- **Cell line stratification**: Individual PDO model responses
- **Concentration ranges**: Full factorial design analysis
- **Statistical validation**: Multiple replicate support

## File Structure
```
pdo_synergy/
├── SynergyFinder.R          # Main analysis script
├── datasets/
│   └── FLV_TBF_data.csv    # Input data file
├── reports/                 # Output directory (auto-generated)
└── README.md               # This file
```

## Notes
- **Cell line exclusion**: Block ID 64 is excluded from analysis
- **Data imputation**: Missing values are handled automatically
- **Output organization**: Results saved in date-stamped folders
- **Python integration**: Some advanced analyses use Python backend

## Dependencies
- R 4.0.0+
- SynergyFinder package
- readr, reticulate packages
- Python environment with synergy analysis tools

## Citation
If using this analysis pipeline, please cite:
- SynergyFinder R package
- Original data sources
- Methodology references
