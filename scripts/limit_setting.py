"""
Complete Limit Setting Analysis - Following Luis's Thesis
Reproduces the methodology from MICROBOONE-NOTE-1118-PUB and Luis's PhD thesis

This script:
1. Loads signal and background CNN scores
2. Creates histograms (like Luis's 20 bins)
3. Builds a pyhf statistical model with uncertainties
4. Calculates CLs limits at 90% confidence level
5. Produces plots and tables like Luis's thesis
"""
# Import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

import pyhf
import json

from scipy import interpolate

## Step 1: Load Data
print("STEP 1: LOADING DATA")

DATA_DIR = Path('/home/paul/Msci/outputs/inference/run1/')

signal_df = pd.read_csv(DATA_DIR / 'run1_dt_ratio_0.6_ma_0.05_pi0_larcv_cropped_scores.csv')
nu_overlay_df = pd.read_csv(DATA_DIR / 'run1_NuMI_nu_overlay_larcv_cropped_scores.csv')
dirt_df = pd.read_csv(DATA_DIR / 'run1_NuMI_dirt_larcv_cropped_scores.csv')
offbeam_df = pd.read_csv(DATA_DIR / 'run1_offbeam_larcv_cropped_full_set_scores.csv')

# Extract CNN scores
signal_scores = signal_df['signal_score'].values
nu_scores = nu_overlay_df['signal_score'].values
dirt_scores = dirt_df['signal_score'].values
offbeam_scores = offbeam_df['signal_score'].values

# Filter out invalid scores (must be in [-1, 1])
signal_scores = signal_scores[(signal_scores >= -1) & (signal_scores <= 1)]
nu_scores = nu_scores[(nu_scores >= -1) & (nu_scores <= 1)]
dirt_scores = dirt_scores[(dirt_scores >= -1) & (dirt_scores <= 1)]
offbeam_scores = offbeam_scores[(offbeam_scores >= -1) & (offbeam_scores <= 1)]

print(f"\nSignal score range: [{signal_scores.min():.3f}, {signal_scores.max():.3f}]")
print(f"Background score range: [{nu_scores.min():.3f}, {nu_scores.max():.3f}]")

# Step 2: Create Histograms
# Define binning (20 bins from 0 to 1)
n_bins = 20
bins = np.linspace(0, 1, n_bins + 1)

# Add weights used by Luis
weight_nu = 0.09
weight_dirt = 0.10
weight_offbeam = 0.66


# Create histograms
signal_hist, _ = np.histogram(signal_scores, bins=bins)
nu_hist, _ = np.histogram(nu_scores, bins=bins, weights=np.ones(len(nu_scores))*weight_nu)
dirt_hist, _ = np.histogram(dirt_scores, bins=bins, weights=np.ones(len(dirt_scores))*weight_dirt)
offbeam_hist, _ = np.histogram(offbeam_scores, bins=bins, weights=np.ones(len(offbeam_scores))*weight_offbeam)

# Total background
background_hist = nu_hist + dirt_hist + offbeam_hist

print(f"\nTotal signal events: {signal_hist.sum()}")
print(f"Total background events: {background_hist.sum()}")

# Build Pyhf model
def build_pyhf_model(signal_hist, nu_hist, dirt_hist, offbeam_hist):
    """
    Build pyhf model
    
    Systematic uncertainties:
    - Signal: POT (2%), meson flux (10%), cross-section (30%), detector (5%)
    - Background: POT (2%), flux, GENIE, detector (5%), statistical
    """
    
    # Convert to lists for JSON
    signal_data = signal_hist.tolist()
    nu_data = nu_hist.tolist()
    dirt_data = dirt_hist.tolist()
    offbeam_data = offbeam_hist.tolist()
    
    # Build model specification
    # NEED TO INCLUDE RUN 2 LATER
    spec = {
        "channels": [
            {
                "name": "CNN_score",
                "samples": [
                    {
                        "name": "signal",
                        "data": signal_data,
                        "modifiers": [
                            # Signal strength (parameter of interest)
                            {
                                "name": "mu_sig",
                                "type": "normfactor",
                                "data": None
                            },
                            # POT scaling uncertainty (±2%)
                            {
                                "name": "POT_signal",
                                "type": "normsys",
                                "data": {"hi": 1.02, "lo": 0.98}
                            },
                            # Meson flux uncertainty (±10%)
                            {
                                "name": "meson_flux",
                                "type": "normsys",
                                "data": {"hi": 1.10, "lo": 0.90}
                            },
                            # Cross-section uncertainty (±30%)
                            {
                                "name": "xsec_signal",
                                "type": "normsys",
                                "data": {"hi": 1.30, "lo": 0.70}
                            },
                            # Detector uncertainty (±5%, correlated shape)
                            {
                                "name": "detector_signal",
                                "type": "normsys",
                                "data": {"hi": 1.05, "lo": 0.95}
                            },
                            # MC statistical uncertainty
                            {
                                "name": "staterror_signal",
                                "type": "staterror",
                                "data": np.sqrt(signal_hist).tolist()
                            }
                        ]
                    },
                    {
                        "name": "nu_overlay",
                        "data": nu_data,
                        "modifiers": [
                            # POT scaling
                            {
                                "name": "POT_background",
                                "type": "normsys",
                                "data": {"hi": 1.02, "lo": 0.98}
                            },
                            # Flux uncertainty (uncorrelated, ~10% per bin)
                            {
                                "name": "flux_nu",
                                "type": "normsys",
                                "data": {"hi": 1.10, "lo": 0.90}
                            },
                            # GENIE cross-section uncertainty (~20%)
                            {
                                "name": "genie_nu",
                                "type": "normsys",
                                "data": {"hi": 1.20, "lo": 0.80}
                            },
                            # Detector uncertainty
                            {
                                "name": "detector_background",
                                "type": "normsys",
                                "data": {"hi": 1.05, "lo": 0.95}
                            },
                            # Statistical uncertainty
                            {
                                "name": "staterror_nu",
                                "type": "staterror",
                                "data": np.sqrt(nu_hist).tolist()
                            }
                        ]
                    },
                    {
                        "name": "dirt",
                        "data": dirt_data,
                        "modifiers": [
                            # Normalization uncertainty (~30%)
                            {
                                "name": "dirt_norm",
                                "type": "normsys",
                                "data": {"hi": 1.30, "lo": 0.70}
                            },
                            # Statistical uncertainty
                            {
                                "name": "staterror_dirt",
                                "type": "staterror",
                                "data": np.sqrt(dirt_hist).tolist()
                            }
                        ]
                    },
                    {
                        "name": "beam_off",
                        "data": offbeam_data,
                        "modifiers": [
                            # Statistical uncertainty only (data-driven)
                            {
                                "name": "staterror_offbeam",
                                "type": "staterror",
                                "data": np.sqrt(offbeam_hist).tolist()
                            }
                        ]
                    }
                ]
            }
        ],
        "observations": [
            {
                "name": "CNN_score",
                "data": background_hist.tolist()  # For expected limit, use background
            }
        ],
        "measurements": [
            {
                "name": "dark_trident_search",
                "config": {
                    "poi": "mu_sig",
                    "parameters": []
                }
            }
        ],
        "version": "1.0.0"
    }
    
    return pyhf.Workspace(spec)

# Build the model
workspace = build_pyhf_model(signal_hist, nu_hist, dirt_hist, offbeam_hist)
model = workspace.model()

print(f"Model built successfully!")
print(f"  - Number of bins: {len(signal_hist)}")
print(f"  - Number of nuisance parameters: {model.config.npars - 1}")
print(f"  - Parameter of interest: {model.config.poi_name}")

# Calculate observed limit (using background as data) (Asimov analysis)
# For expected limit, we use background as observation (In real analysis, you'd use actual data here)
observations = background_hist.tolist()
data = observations + model.config.auxdata # Adds constraints

print("Scanning signal strength μ values...")

# Scan over mu values
mu_scan = np.linspace(0.001, 3.0, 50)
CLs_obs_values = []
CLs_exp_values = []

for i, mu_test in enumerate(mu_scan):
    if i % 10 == 0:
        print(f"  Testing μ = {mu_test:.3f}...") # Print only 10th iteration
    # Calculate CLs using asymptotic approximation (the bread and butter)
    CLs_obs, CLs_exp_band = pyhf.infer.hypotest(
        mu_test,
        data,
        model,
        test_stat="qtilde",
        return_expected_set=True,
        calctype="asymptotics"  
    )
    
    CLs_obs_values.append(CLs_obs)
    CLs_exp_values.append(CLs_exp_band[2])  # Median expected, Same as observed for Asimov

CLs_obs_values = np.array(CLs_obs_values)
CLs_exp_values = np.array(CLs_exp_values)

# Find where CLs crosses 0.10 (90% CL)
def find_limit(mu_values, CLs_values, alpha=0.10):
    """Find mu where CLs = alpha by interpolation"""
    # Reverse arrays for interpolation (CLs decreases with mu)
    if CLs_values[-1] > CLs_values[0]:
        # Already in right order
        pass
    else:
        mu_values = mu_values[::-1]
        CLs_values = CLs_values[::-1]
    
    # Find crossing point
    if CLs_values[0] < alpha:
        return mu_values[0]  # Limit below scan range
    if CLs_values[-1] > alpha:
        return mu_values[-1]  # Limit above scan range
    
    # Interpolate
    limit = np.interp(alpha, CLs_values[::-1], mu_values[::-1])
    return limit

mu_limit_obs = find_limit(mu_scan, CLs_obs_values, alpha=0.10)
mu_limit_exp = find_limit(mu_scan, CLs_exp_values, alpha=0.10)

print(f"\n90% CL Upper Limit on Signal Strength μ:")
print(f"  Expected limit: μ < {mu_limit_exp:.3f}")
print(f"  Observed limit: μ < {mu_limit_obs:.3f}")

# Translate to number of events
n_signal_total = signal_hist.sum()
n_signal_excluded = mu_limit_exp * n_signal_total

print(f"\nExcluded number of signal events:")
print(f"  {n_signal_excluded:.1f} events at 90% CL")
print(f"  (Benchmark signal has {n_signal_total:.0f} events)")

# Translate to epsilon^2
benchmark_eps2 = 1e-6  # Benchmark coupling from simulation
M_A_prime = 0.05  # GeV (50 MeV)
alpha_D = 0.1

eps2_limit_exp = benchmark_eps2 * mu_limit_exp
eps2_limit_obs = benchmark_eps2 * mu_limit_obs

print(f"\nPhysics interpretation:")
print(f"  For M_A' = {M_A_prime*1000:.0f} MeV, α_D = {alpha_D}:")
print(f"  Expected: ε² < {eps2_limit_exp:.2e} at 90% CL")
print(f"  Observed: ε² < {eps2_limit_obs:.2e} at 90% CL")

# Compare to Luis's Table 11.5 value for 50 MeV
luis_limit = 1.19e-7
print(f"\n  Luis's expected: ε² < {luis_limit:.2e}")
print(f"  Calculated expected limit: ε² < {eps2_limit_exp:.2e}")
print(f"  Ratio: {eps2_limit_exp/luis_limit:.2f}x")

# Step 5: Making plots

# Plot 1: CLs vs μ (Luis's Figure 6.4)
fig, ax = plt.subplots(figsize=(10, 7))

ax.plot(mu_scan, CLs_obs_values, 'b-', linewidth=2, label='$CL_s$ (observed)')
ax.plot(mu_scan, CLs_exp_values, 'b--', linewidth=2, label='$CL_s$ (expected)')
ax.axhline(0.10, color='r', linestyle='-', linewidth=2, label='α = 0.10 (90% CL)')
ax.axvline(mu_limit_exp, color='g', linestyle=':', alpha=0.7, label=f'Expected limit = {mu_limit_exp:.2f}')
ax.axvline(mu_limit_obs, color='b', linestyle=':', alpha=0.7, label=f'Observed limit = {mu_limit_obs:.2f}')

ax.set_xlabel('Signal Strength μ', fontsize=14)
ax.set_ylabel('$CL_s$', fontsize=14)
ax.set_title('CLs vs Signal strength μ (90% CL Limit)', fontsize=16)
ax.set_ylim(0, 1.0)
ax.set_xlim(0, mu_scan[-1])
ax.grid(True, alpha=0.3)
ax.legend(fontsize=11, loc='upper right')

plt.tight_layout()
plt.savefig('/home/paul/Msci/outputs/limit_setting')
print("Saved cls_limit_scan.png")
plt.close()

# Step 6: Final results
results_table = (
    "90% CL Limits on Dark Trident Signal\n"
    "(Following Luis's Thesis Methodology)\n"
    "\n"
    f"Dark Photon Mass (M_A'):      {M_A_prime*1000:.0f} MeV\n"
    f"Dark Coupling (alpha_D):      {alpha_D}\n"
    "Mass Ratio (M_chi / M_A'):    0.6\n"
    "\n"
    "SIGNAL STRENGTH LIMITS\n"
    f"  Expected limit:   mu < {mu_limit_exp:.3f}\n"
    f"  Observed limit:   mu < {mu_limit_obs:.3f}\n"
    "\n"
    "EXCLUDED SIGNAL EVENTS\n"
    f"  Expected excluded events:   > {n_signal_excluded:.1f}\n"
    f"  Benchmark signal events:      {n_signal_total:.0f}\n"
    f"  Signal excluded?              {n_signal_total > n_signal_excluded}\n"
    "\n"
    "PHYSICS PARAMETERS (eps^2)\n"
    f"  Expected limit:   eps^2 < {eps2_limit_exp:.2e}\n"
    f"  Observed limit:   eps^2 < {eps2_limit_obs:.2e}\n"
    f"  Luis's limit:     eps^2 < {luis_limit:.2e}\n"
    f"  Ratio (you/Luis): {eps2_limit_exp/luis_limit:.2f}x\n"
    "\n"
    "STATISTICS SUMMARY\n"
    f"  Background total: {background_hist.sum():.0f}\n"
    f"    - Nu overlay:   {nu_hist.sum():.0f}\n"
    f"    - Dirt:         {dirt_hist.sum():.0f}\n"
    f"    - Beam-off:     {offbeam_hist.sum():.0f}\n"
    "\n"
    f"  Signal events:    {signal_hist.sum():.0f}\n"
    f"  S/B ratio:        {signal_hist.sum()/background_hist.sum():.3f}\n"
)

print(results_table)

