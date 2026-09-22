# FLV-TBF Dataset Overview

## Dataset Summary
- **Total Measurements**: 2,150 data points
- **Primary Drug**: Fluvastatin (FLV) - HMG-CoA reductase inhibitor (statin)
- **Secondary Drug**: Terbinafine (TBF) - Antifungal agent
- **Cell Lines**: 7 different breast cancer patient-derived organoid (PDO) models
- **Experimental Repeats**: Up to 5 replicates per condition
- **Drug Combination**: FLV + TBF (statin + antifungal combination therapy)

## Cell Line Coverage

| Block ID | Measurements | Repeats | FLV Conc. | TBF Conc. | Mean Response (%) | SD Response (%) |
|----------|--------------|---------|------------|-----------|------------------|-----------------|
| 58       | 340          | 3       | 18         | 18        | 45.2             | 28.4            |
| 66       | 180          | 2       | 18         | 18        | 42.8             | 29.1            |
| 95       | 540          | 5       | 18         | 18        | 38.9             | 31.2            |
| 132      | 240          | 2       | 18         | 18        | 41.3             | 30.7            |
| 137      | 240          | 2       | 18         | 18        | 39.7             | 32.1            |
| 143      | 240          | 2       | 18         | 18        | 40.2             | 29.8            |
| 180      | 250          | 3       | 18         | 18        | 43.1             | 30.5            |

## Drug Concentration Ranges

### Fluvastatin (FLV) Concentrations
**18 different concentration levels tested:**
- **Low range**: 0, 0.1, 0.17, 0.3, 0.5, 0.9 μM
- **Medium range**: 1.6, 2.2, 2.7, 3.3, 4.8, 5, 7.5, 8.3 μM
- **High range**: 11.1, 14.4, 16.7, 25 μM

### Terbinafine (TBF) Concentrations
**18 different concentration levels tested:**
- **Low range**: 0, 0.1, 0.2, 0.4, 0.5, 0.83, 1, 1.6, 2.2 μM
- **Medium range**: 3.15, 4.9, 6.3, 8.9, 10.6, 12.5 μM
- **High range**: 17.8, 25, 50 μM

## Experimental Design

### Concentration Matrix
- **FLV concentrations**: 18 levels (0-25 μM)
- **TBF concentrations**: 18 levels (0-50 μM)
- **Total combinations**: 324 concentration pairs
- **Experimental format**: Full factorial design matrix
- **Control conditions**: Single drug (FLV only, TBF only) and combination treatments

### Quality Control
- **Replicate coverage**: 2-5 experimental repeats per condition
- **Cell lines**: 7 different PDO models for biological validation
- **Concentration units**: Micromolar (μM) for both drugs
- **Response measurement**: Cell viability assays

## Data Structure

### Variables
- **`date`**: Experiment date (Excel date format)
- **`conc1`**: FLV concentration (μM)
- **`conc2`**: TBF concentration (μM)
- **`response`**: Cell viability measurement (%)
- **`repeats`**: Experimental replicate number (1-5)
- **`drug1`**: Primary drug identifier (FLV)
- **`drug2`**: Secondary drug identifier (TBF)
- **`conc_unit1`**: Concentration unit for FLV (μM)
- **`conc_unit2`**: Concentration unit for TBF (μM)
- **`block_id`**: Cell line/patient sample identifier

## Notes
- Dataset represents high-throughput screening data from combination therapy studies
- All measurements are from standardized cell viability assays
- Concentration ranges cover therapeutic to toxic levels
- Data suitable for advanced statistical analysis, modeling, and machine learning
- Results can inform clinical trial design and personalized medicine approaches
