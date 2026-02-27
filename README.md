# Blood-Stage Minibinder EC50 Analysis

EC50 calculation pipeline for anti-PfRH5 minibinder blood-stage inhibition screen data.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run analysis
python calculate_ec50.py
```

## Input Data

Two Excel files with biological replicates:
- `BloodstageMinibinderScreen_12-12-25.xlsx` (Rep 1)
- `BloodstageMinibinderScreen_Rep2_12-19-25.xlsx` (Rep 2)

Each file contains 384-well plate data with luminescence readings from a Dd2-Luc blood-stage growth inhibition assay.

## Output Files

- `EC50_bar_graph.png` - Publication-quality bar graph (300 dpi)
- `EC50_bar_graph.pdf` - PDF version
- `EC50_values_table.csv` - EC50 values for all compounds

## EC50 Calculation Method

### 1. Data Extraction
Raw luminescence data is extracted from each Excel file. Each compound has 11 concentration points with duplicate wells (22 data points per compound per replicate).

### 2. Normalization
PBS control wells serve as the "max signal" baseline (100% parasite growth, 0% inhibition):

```
% inhibition = 100 × (1 - signal / PBS_mean)
```

### 3. Curve Fitting
A 4-parameter logistic (4PL) model is fit to each dose-response curve:

```
y = bottom + (top - bottom) / (1 + (EC50/x)^hill)
```

This is the same model used by GraphPad Prism. Parameters:
- `bottom`: minimum response (% inhibition at low concentration)
- `top`: maximum response (% inhibition at high concentration)
- `EC50`: concentration producing 50% effect
- `hill`: Hill slope (steepness of curve)

### 4. Unit Conversion

| Compound Type | Original Unit | Conversion |
|---------------|---------------|------------|
| Drug controls (ATQ, GNF179) | µM | × 1000 → nM |
| Minibinders (~16 kDa) | mg/mL | × 10⁹ / 16000 → nM |
| Monoclonal Abs (~150 kDa) | mg/mL | × 10⁹ / 150000 → nM |

### 5. Replicate Handling
- EC50 calculated separately for each biological replicate
- Geometric mean used for final EC50 (standard for log-distributed values)
- Individual replicate points shown on bar graph

## Compound Categories

| Category | Compounds | Color |
|----------|-----------|-------|
| Drug Control | ATQ, GNF179 | Red |
| Monoclonal Ab | R5.004, Cy.003 | Green |
| Minibinder | A8, B7, B9, B8, B4, B10, B11, B6, B1 | Blue |

## Data Quality Notes

- **A4**: Excluded - did not achieve >30% inhibition at highest concentration
- **B6 Rep 1**: Excluded - insufficient inhibition
- **B7 Rep 2**: Wells H13-23 excluded due to Integra pipetting error (noted in original file)

## Results Summary

### Most Potent Minibinder
**A8**: EC50 = 79 nM

### Potent Minibinders (100-250 nM)
- Cy.003: 208 nM (mAb control)
- B7: 210 nM
- B9: 236 nM

### Drug Controls (reference)
- GNF179: 0.6 nM
- ATQ: 0.9 nM

## Dependencies

- Python 3.8+
- pandas
- numpy
- scipy
- matplotlib
- openpyxl (for Excel reading)

## Usage

```bash
# Default (uses files in current directory)
python calculate_ec50.py

# Specify files
python calculate_ec50.py --rep1 path/to/rep1.xlsx --rep2 path/to/rep2.xlsx

# Specify output directory
python calculate_ec50.py --output results/
```
