#!/usr/bin/env python3
"""
EC50 Calculation for Anti-PfRH5 Minibinder Blood-Stage Screen

This script calculates EC50 values from dose-response data in the blood-stage
malaria inhibition assay. It processes Excel files containing luminescence
readings and fits 4-parameter logistic curves to determine EC50.

Author: Jinich Lab
Date: December 2025
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Patch
import warnings
import argparse
import os

warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

# Molecular weights (Da)
MW_MINIBINDER = 16000    # ~16 kDa
MW_MONOCLONAL = 150000   # ~150 kDa

# Compound classifications
DRUG_CONTROLS = ['ATQ', 'GNF179']
MONOCLONAL_ABS = ['R5.004', 'Cy.003']
# All others are assumed to be minibinders

# Known problematic wells to exclude (from experimental notes)
BAD_WELLS_REP2 = [
    'H13', 'H14', 'H15', 'H16', 'H17', 'H18', 'H19', 'H20', 'H21', 'H22', 'H23',
    'I9', 'I10', 'I11'
]

# Plot colors
COLORS = {
    'Drug Control': '#E74C3C',    # Red
    'Minibinder': '#3498DB',      # Blue
    'Monoclonal Ab': '#27AE60',   # Green
}


# =============================================================================
# DATA EXTRACTION FUNCTIONS
# =============================================================================

def extract_data(df):
    """
    Extract compound-concentration-signal data from Excel format.

    The Excel files have a specific layout where columns 26-31 contain:
    - Column 26: Plate Name
    - Column 27: Well Location
    - Column 28: Compound name
    - Column 29: Batch Name
    - Column 30: Concentration
    - Column 31: Signal (luminescence)

    Parameters
    ----------
    df : pandas.DataFrame
        Raw Excel data read with header=None

    Returns
    -------
    pandas.DataFrame
        Cleaned data with columns: compound, concentration, signal, well
    """
    data = []
    for i in range(8, len(df)):
        row = df.iloc[i]
        cmpd = row[28]
        conc = row[30]
        signal = row[31]
        well = row[27]

        # Skip empty/control entries
        if pd.notna(cmpd) and pd.notna(conc) and pd.notna(signal):
            if cmpd not in ['CMPD', 'Empty', 'PBS']:
                try:
                    data.append({
                        'compound': str(cmpd),
                        'concentration': float(conc),
                        'signal': float(signal),
                        'well': str(well) if pd.notna(well) else ''
                    })
                except (ValueError, TypeError):
                    pass

    return pd.DataFrame(data)


def get_pbs_mean(df):
    """
    Get mean signal from PBS control wells (no inhibition baseline).

    Parameters
    ----------
    df : pandas.DataFrame
        Raw Excel data

    Returns
    -------
    float
        Mean luminescence signal from PBS wells
    """
    signals = []
    for i in range(8, len(df)):
        row = df.iloc[i]
        cmpd = row[28]
        signal = row[31]
        if cmpd == 'PBS' and pd.notna(signal):
            signals.append(float(signal))
    return np.mean(signals) if signals else 0


# =============================================================================
# EC50 CALCULATION
# =============================================================================

def four_param_logistic(x, bottom, top, ec50, hill):
    """
    4-parameter logistic function for dose-response curves.

    y = bottom + (top - bottom) / (1 + (x/ec50)^hill)

    This is the standard model used by GraphPad Prism for dose-response fitting.
    """
    return bottom + (top - bottom) / (1 + (x / ec50) ** hill)


def calculate_ec50(conc, sig, max_signal, compound_name, min_inhibition=30):
    """
    Calculate EC50 from dose-response data using 4-parameter logistic fit.

    Parameters
    ----------
    conc : array-like
        Concentration values
    sig : array-like
        Signal (luminescence) values
    max_signal : float
        Maximum signal (from PBS controls) for normalization
    compound_name : str
        Name of compound (for error messages)
    min_inhibition : float
        Minimum % inhibition required to attempt EC50 calculation

    Returns
    -------
    tuple
        (ec50_value, status_message)
    """
    if len(conc) < 4:
        return np.nan, "Insufficient data points"

    # Sort by concentration
    sort_idx = np.argsort(conc)
    conc = np.array(conc)[sort_idx]
    sig = np.array(sig)[sort_idx]

    # Convert to % inhibition
    # Higher signal = more parasite growth = less inhibition
    pct_inhibition = 100 * (1 - sig / max_signal)
    max_inhib = np.max(pct_inhibition)

    # Check if there's meaningful inhibition
    if max_inhib < min_inhibition:
        return np.nan, f"Max inhibition only {max_inhib:.1f}%"

    try:
        # Fit inhibition curve
        # Using form: y = bottom + (top-bottom) / (1 + (EC50/x)^hill)
        def inhib_curve(x, ec50, hill, bottom, top):
            return bottom + (top - bottom) / (1 + (ec50 / x) ** hill)

        popt, pcov = curve_fit(
            inhib_curve,
            conc,
            pct_inhibition,
            p0=[np.median(conc), 1, 0, 100],
            bounds=(
                [conc.min()/1000, 0.1, -20, 50],
                [conc.max()*1000, 10, 50, 120]
            ),
            maxfev=10000
        )

        ec50 = popt[0]

        # Validate EC50 is within reasonable range
        if ec50 < conc.min() / 10:
            return ec50, "EC50 below tested range"
        if ec50 > conc.max() * 10:
            return np.nan, "EC50 above tested range"

        return ec50, "OK"

    except Exception as e:
        return np.nan, f"Fit failed: {str(e)}"


def assign_category(compound):
    """Assign compound to category based on predefined lists."""
    if compound in DRUG_CONTROLS:
        return 'Drug Control'
    elif compound in MONOCLONAL_ABS:
        return 'Monoclonal Ab'
    else:
        return 'Minibinder'


def convert_to_nM(ec50_value, category):
    """
    Convert EC50 from original units to nM.

    Drug controls: µM → nM (×1000)
    Minibinders: mg/mL → nM using MW = 16 kDa
    Monoclonal Abs: mg/mL → nM using MW = 150 kDa
    """
    if pd.isna(ec50_value):
        return np.nan

    if category == 'Drug Control':
        # µM to nM
        return ec50_value * 1000
    elif category == 'Monoclonal Ab':
        # mg/mL to nM: nM = (mg/mL × 1e9) / MW
        return (ec50_value * 1e9) / MW_MONOCLONAL
    else:
        # Minibinder: mg/mL to nM
        return (ec50_value * 1e9) / MW_MINIBINDER


# =============================================================================
# PLOTTING
# =============================================================================

def create_bar_plot(results_df, output_prefix='EC50_bar_graph'):
    """
    Create publication-quality bar plot of EC50 values.

    Parameters
    ----------
    results_df : pandas.DataFrame
        DataFrame with EC50 results
    output_prefix : str
        Prefix for output files (will create .png and .pdf)
    """
    # Set publication-quality defaults
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2

    # Filter and sort
    plot_df = results_df[results_df['ec50_mean_nM'].notna()].copy()
    plot_df = plot_df.sort_values('ec50_mean_nM').reset_index(drop=True)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Bar positions
    x = np.arange(len(plot_df))
    bar_width = 0.7

    # Plot bars
    bar_colors = [COLORS[cat] for cat in plot_df['category']]
    bars = ax.bar(x, plot_df['ec50_mean_nM'], width=bar_width, color=bar_colors,
                  edgecolor='black', linewidth=1, alpha=0.8, zorder=2)

    # Overlay individual replicate data points
    for i, (_, row) in enumerate(plot_df.iterrows()):
        points = []
        if pd.notna(row['ec50_rep1_nM']):
            points.append(row['ec50_rep1_nM'])
        if pd.notna(row['ec50_rep2_nM']):
            points.append(row['ec50_rep2_nM'])

        jitter = np.linspace(-0.12, 0.12, len(points))
        for pt, jit in zip(points, jitter):
            ax.scatter(i + jit, pt, s=40, c='black', marker='o', zorder=3,
                      edgecolors='white', linewidth=0.5)

    # Add reference line at 100 nM
    ax.axhline(y=100, color='gray', linestyle='--', linewidth=1.5, zorder=1, alpha=0.7)
    ax.text(len(plot_df) - 0.5, 120, '100 nM', fontsize=9, color='gray', ha='right', va='bottom')

    # Set log scale
    ax.set_yscale('log')

    # Axis settings
    ax.set_xlim(-0.6, len(plot_df) - 0.4)
    ax.set_ylim(0.1, 200000)
    ax.set_xlabel('Compound', fontsize=12, fontweight='bold')
    ax.set_ylabel('EC$_{50}$ (nM)', fontsize=12, fontweight='bold')

    # X-tick labels
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df['compound'], rotation=45, ha='right', fontsize=10)

    # Clean up spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)

    # Legend
    legend_elements = [
        Patch(facecolor=COLORS['Drug Control'], edgecolor='black', label='Drug Control'),
        Patch(facecolor=COLORS['Monoclonal Ab'], edgecolor='black', label='Monoclonal Ab'),
        Patch(facecolor=COLORS['Minibinder'], edgecolor='black', label='Minibinder'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=False, fontsize=10)

    plt.tight_layout()

    # Save
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_prefix}.pdf', bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved: {output_prefix}.png (300 dpi)")
    print(f"Saved: {output_prefix}.pdf")


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def analyze_screen(file1, file2, output_dir='.'):
    """
    Run full EC50 analysis on screening data.

    Parameters
    ----------
    file1 : str
        Path to first replicate Excel file
    file2 : str
        Path to second replicate Excel file
    output_dir : str
        Directory for output files

    Returns
    -------
    pandas.DataFrame
        Results table with EC50 values
    """
    print("=" * 70)
    print("EC50 ANALYSIS FOR BLOOD-STAGE MINIBINDER SCREEN")
    print("=" * 70)

    # Read data
    print(f"\nReading {file1}...")
    df1 = pd.read_excel(file1, sheet_name='Sheet1', header=None)
    print(f"Reading {file2}...")
    df2 = pd.read_excel(file2, sheet_name='Sheet1', header=None)

    # Extract dose-response data
    data1 = extract_data(df1)
    data2 = extract_data(df2)

    # Clean Rep2 data (remove known bad wells)
    data2_clean = data2[~data2['well'].isin(BAD_WELLS_REP2)]

    # Get max signals for normalization
    max_sig1 = get_pbs_mean(df1)
    max_sig2 = get_pbs_mean(df2)

    print(f"\nPBS control signals: Rep1={max_sig1:.0f}, Rep2={max_sig2:.0f}")

    # Get all compounds
    compounds = sorted(set(data1['compound'].unique()) | set(data2['compound'].unique()))
    compounds = [c for c in compounds if c not in ['CMPD', 'Empty', 'PBS']]

    print(f"\nAnalyzing {len(compounds)} compounds...")

    # Calculate EC50 for each compound
    results = []

    for compound in compounds:
        # Rep 1
        cmp_data1 = data1[data1['compound'] == compound]
        if len(cmp_data1) > 0:
            ec50_1, status_1 = calculate_ec50(
                cmp_data1['concentration'].values,
                cmp_data1['signal'].values,
                max_sig1, compound
            )
        else:
            ec50_1, status_1 = np.nan, "No data"

        # Rep 2 (cleaned)
        cmp_data2 = data2_clean[data2_clean['compound'] == compound]
        if len(cmp_data2) > 0:
            ec50_2, status_2 = calculate_ec50(
                cmp_data2['concentration'].values,
                cmp_data2['signal'].values,
                max_sig2, compound
            )
        else:
            ec50_2, status_2 = np.nan, "No data"

        category = assign_category(compound)

        results.append({
            'compound': compound,
            'category': category,
            'ec50_rep1': ec50_1,
            'ec50_rep2': ec50_2,
            'status_rep1': status_1,
            'status_rep2': status_2
        })

    results_df = pd.DataFrame(results)

    # Convert to nM
    for idx, row in results_df.iterrows():
        ec50_1_nM = convert_to_nM(row['ec50_rep1'], row['category'])
        ec50_2_nM = convert_to_nM(row['ec50_rep2'], row['category'])

        results_df.at[idx, 'ec50_rep1_nM'] = ec50_1_nM
        results_df.at[idx, 'ec50_rep2_nM'] = ec50_2_nM

        # Geometric mean of valid replicates
        valid_vals = [v for v in [ec50_1_nM, ec50_2_nM] if pd.notna(v)]
        if len(valid_vals) == 2:
            results_df.at[idx, 'ec50_mean_nM'] = np.sqrt(valid_vals[0] * valid_vals[1])
        elif len(valid_vals) == 1:
            results_df.at[idx, 'ec50_mean_nM'] = valid_vals[0]
        else:
            results_df.at[idx, 'ec50_mean_nM'] = np.nan

    # Print results
    print("\n" + "=" * 70)
    print("EC50 VALUES (nM)")
    print("=" * 70)
    print(f"\n{'Compound':<10} {'Category':<14} {'Rep1':>12} {'Rep2':>12} {'Mean':>12}")
    print("-" * 62)

    for _, row in results_df.sort_values('ec50_mean_nM').iterrows():
        ec50_1 = f"{row['ec50_rep1_nM']:.1f}" if pd.notna(row['ec50_rep1_nM']) else "N/A"
        ec50_2 = f"{row['ec50_rep2_nM']:.1f}" if pd.notna(row['ec50_rep2_nM']) else "N/A"
        mean = f"{row['ec50_mean_nM']:.1f}" if pd.notna(row['ec50_mean_nM']) else "N/A"
        print(f"{row['compound']:<10} {row['category']:<14} {ec50_1:>12} {ec50_2:>12} {mean:>12}")

    # Save results
    output_csv = os.path.join(output_dir, 'EC50_values_table.csv')
    results_df.to_csv(output_csv, index=False)
    print(f"\nSaved: {output_csv}")

    # Create plot
    create_bar_plot(results_df, os.path.join(output_dir, 'EC50_bar_graph'))

    return results_df


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculate EC50 values from blood-stage malaria screen data'
    )
    parser.add_argument(
        '--rep1',
        default='BloodstageMinibinderScreen_12-12-25.xlsx',
        help='Path to replicate 1 Excel file'
    )
    parser.add_argument(
        '--rep2',
        default='BloodstageMinibinderScreen_Rep2_12-19-25.xlsx',
        help='Path to replicate 2 Excel file'
    )
    parser.add_argument(
        '--output', '-o',
        default='.',
        help='Output directory for results'
    )

    args = parser.parse_args()

    results = analyze_screen(args.rep1, args.rep2, args.output)
