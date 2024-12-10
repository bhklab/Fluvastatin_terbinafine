# Fluvastatin_terbinafine Combo Analysis
 
## Introduction
This is a project to analyze combo data in this manuscript:

Genome-wide CRISPR Screen Identifies Squalene Epoxidase as a Fluvastatin Sensitizer to Effectively Target the Mevalonate Pathway in Breast Cancer, 2024, Leeuwen et. al..

The script for the analysis is written in R.


----

## Dependencies:
These R packages need to be installed in order to run the analysis script:
- PharmacoGx
- magicaxis
- abind
- robustbase
- Biobase
- synergyfinder (1.8.0)
- ComplexHeatmap
- circlize
- ggplot2
- ggpubr
- reshape
- snowfall
- GSA
- piano
- scales
- ggrepel
- dplyr
- tidyr
- BoutrosLab.plotting.general
- UpSetR


----
## Reproducibility of the Analysis:
- Inside the main directory, there is an R script file named "combo_analysis.R". Running this script will regenerate the combo data and analyses in the manuscript linked above.

**Important Note:** the user needs to set the working directory inside the script file before running it, i.e. setting the working directory to "combo_analysis"



