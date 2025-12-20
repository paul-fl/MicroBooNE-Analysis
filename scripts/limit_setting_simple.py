"""
Simple Limit Setting with pyhf
------------------------------
Minimal working example - no fancy systematics, just the basics.

What this does:
1. Load CNN scores from CSV files
2. Build histograms
3. Create a simple pyhf model (signal + background)
4. Calculate 90% CL upper limit on signal strength mu
5. Plot CLs vs mu (Brazil band)
"""

import numpy as np
import pandas as pd
import pyhf
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION - Edit these paths for your setup
# =============================================================================

# File paths - UPDATE THESE TO YOUR ACTUAL PATHS
RUN1_DIR = "/home/paul/Msci/outputs/inference/run1/"
RUN3_DIR = "/home/paul/Msci/outputs/inference/run3/"

# Files to load
FILES = {
    "run1": {
        "signal": "run1_dt_ratio_0.6_ma_0.05_pi0_larcv_cropped_scores.csv",
        "nu": "run1_NuMI_nu_overlay_larcv_cropped_scores.csv",
        "dirt": "run1_NuMI_dirt_larcv_cropped_scores.csv",
        "offbeam": "run1_offbeam_larcv_cropped_full_set_scores.csv",
    },
    "run3": {
        "signal": "run3_dt_ratio_0.6_ma_0.05_pi0_larcv_cropped_scores.csv",
        "nu": "run3_nu_overlay_larcv_cropped_scores.csv",
        "dirt": "run3_dirt_larcv_cropped_scores.csv",
        "offbeam": "run3_offbeam_larcv_cropped_scores.csv",
    },
}

# Weights from Luis's Table 7.4
WEIGHTS = {
    "run1": {"nu": 0.09, "dirt": 0.10, "offbeam": 0.66},
    "run3": {"nu": 0.24, "dirt": 0.17, "offbeam": 0.30},
}

# =============================================================================
# STEP 1: Load CSV files and extract scores
# =============================================================================

def load_scores(directory, filename):
    """Load CNN signal scores from a CSV file."""
    df = pd.read_csv(directory + filename)
    scores = df["signal_score"].values
    # Filter to valid range [0, 1]
    scores = scores[(scores >= 0) & (scores <= 1)]
    return scores


print("Loading data...")

# Load Run 1
sig1 = load_scores(RUN1_DIR, FILES["run1"]["signal"])
nu1 = load_scores(RUN1_DIR, FILES["run1"]["nu"])
dirt1 = load_scores(RUN1_DIR, FILES["run1"]["dirt"])
off1 = load_scores(RUN1_DIR, FILES["run1"]["offbeam"])

# Load Run 3
sig3 = load_scores(RUN3_DIR, FILES["run3"]["signal"])
nu3 = load_scores(RUN3_DIR, FILES["run3"]["nu"])
dirt3 = load_scores(RUN3_DIR, FILES["run3"]["dirt"])
off3 = load_scores(RUN3_DIR, FILES["run3"]["offbeam"])

print(f"  Run1 - Signal: {len(sig1)}, Nu: {len(nu1)}, Dirt: {len(dirt1)}, Offbeam: {len(off1)}")
print(f"  Run3 - Signal: {len(sig3)}, Nu: {len(nu3)}, Dirt: {len(dirt3)}, Offbeam: {len(off3)}")

# =============================================================================
# STEP 2: Build histograms
# =============================================================================

print("\nBuilding histograms...")

# Binning: 20 bins from 0 to 1
bins = np.linspace(0, 1, 21)

def make_hist(scores, weight=1.0):
    """Create a weighted histogram."""
    counts, _ = np.histogram(scores, bins=bins)
    return counts.astype(float) * weight


# Run 1 histograms (with weights)
h_sig1 = make_hist(sig1, weight=1.0)  # Signal unweighted for now
h_nu1 = make_hist(nu1, weight=WEIGHTS["run1"]["nu"])
h_dirt1 = make_hist(dirt1, weight=WEIGHTS["run1"]["dirt"])
h_off1 = make_hist(off1, weight=WEIGHTS["run1"]["offbeam"])
h_bkg1 = h_nu1 + h_dirt1 + h_off1

# Run 3 histograms (with weights)
h_sig3 = make_hist(sig3, weight=1.0)
h_nu3 = make_hist(nu3, weight=WEIGHTS["run3"]["nu"])
h_dirt3 = make_hist(dirt3, weight=WEIGHTS["run3"]["dirt"])
h_off3 = make_hist(off3, weight=WEIGHTS["run3"]["offbeam"])
h_bkg3 = h_nu3 + h_dirt3 + h_off3

# Combine Run1 + Run3 (concatenate bins)
signal_hist = np.concatenate([h_sig1, h_sig3])
background_hist = np.concatenate([h_bkg1, h_bkg3])

# For Asimov (expected limit): data = background
data_hist = background_hist.copy()

print(f"  Total signal events: {signal_hist.sum():.1f}")
print(f"  Total background events: {background_hist.sum():.1f}")
print(f"  S/B ratio: {signal_hist.sum() / background_hist.sum():.3f}")

# =============================================================================
# STEP 3: Build the pyhf model (SIMPLE VERSION)
# =============================================================================

print("\nBuilding pyhf model...")

# Use pyhf's simple model helper - this is the easiest way!
# It creates a model with:
#   - One signal sample (scaled by mu)
#   - One background sample (with uncorrelated bin-by-bin uncertainties)

# Background uncertainty: just use sqrt(N) for now (Poisson-like)
bkg_uncertainty = np.sqrt(background_hist)

# Handle zero bins (avoid division issues)
bkg_uncertainty = np.where(bkg_uncertainty > 0, bkg_uncertainty, 0.1)

model = pyhf.simplemodels.uncorrelated_background(
    signal=signal_hist.tolist(),
    bkg=background_hist.tolist(),
    bkg_uncertainty=bkg_uncertainty.tolist(),
)

# Build the observation (data + auxiliary data for nuisance parameters)
observations = data_hist.tolist() + model.config.auxdata

print(f"  Model has {model.config.npars} parameters")
print(f"  Parameter of interest: {model.config.poi_name}")

# =============================================================================
# STEP 4: Calculate the upper limit
# =============================================================================

print("\nCalculating upper limit...")

# Scan mu values from 0 to 10
poi_values = np.linspace(0.1, 10.0, 50)

# Use the NEW pyhf API (not the deprecated one)
obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(
    observations,
    model,
    poi_values,
    level=0.10,  # 90% CL (1 - 0.90 = 0.10)
    return_results=True,
)

print("\n" + "=" * 50)
print("RESULTS (90% CL Upper Limits)")
print("=" * 50)
print(f"  Observed limit:    mu < {obs_limit:.3f}")
print(f"  Expected -2sigma:  mu < {exp_limits[0]:.3f}")
print(f"  Expected -1sigma:  mu < {exp_limits[1]:.3f}")
print(f"  Expected (median): mu < {exp_limits[2]:.3f}")
print(f"  Expected +1sigma:  mu < {exp_limits[3]:.3f}")
print(f"  Expected +2sigma:  mu < {exp_limits[4]:.3f}")

# Convert to epsilon^2 (physics parameter)
# Using: eps^2 = (nominal_eps)^2 * sqrt(mu_limit)
nominal_eps = 1e-3
eps2_obs = (nominal_eps**2) * np.sqrt(obs_limit)
eps2_exp = (nominal_eps**2) * np.sqrt(exp_limits[2])

print(f"\n  Observed epsilon^2 < {eps2_obs:.2e}")
print(f"  Expected epsilon^2 < {eps2_exp:.2e}")

# =============================================================================
# STEP 5: Make the Brazil band plot
# =============================================================================

print("\nMaking Brazil band plot...")

from pyhf.contrib.viz import brazil

fig, ax = plt.subplots(figsize=(10, 7))

# Plot the Brazil band
brazil.plot_results(poi_values, results, test_size=0.10, ax=ax)

ax.set_xlabel(r"Signal strength $\mu$", fontsize=14)
ax.set_ylabel(r"$\mathrm{CL}_s$", fontsize=14)
ax.set_title("90% CL Upper Limit (Brazil Band)", fontsize=16)
ax.legend(fontsize=12)

# Set reasonable axis limits
ax.set_xlim(0, min(obs_limit * 3, poi_values[-1]))
ax.set_ylim(0, 1)

plt.tight_layout()
plt.savefig("limit_brazil_band.png", dpi=150)
print("  Saved: limit_brazil_band.png")
plt.show()

# =============================================================================
# STEP 6: Also make a simple CLs scan plot
# =============================================================================

fig, ax = plt.subplots(figsize=(10, 7))

# Extract CLs values from results
# results is a list of (CLs_obs, CLs_exp_band) tuples
cls_obs = np.array([r[0] for r in results])
cls_exp = np.array([r[1][2] for r in results])  # median expected

ax.plot(scan, cls_obs, "k-", linewidth=2, label="Observed CLs")
ax.plot(scan, cls_exp, "k--", linewidth=2, label="Expected CLs")
ax.axhline(0.10, color="r", linestyle="-", linewidth=2, label=r"$\alpha = 0.10$ (90% CL)")
ax.axvline(obs_limit, color="b", linestyle=":", alpha=0.7, label=f"Limit = {obs_limit:.2f}")

ax.set_xlabel(r"Signal strength $\mu$", fontsize=14)
ax.set_ylabel(r"$\mathrm{CL}_s$", fontsize=14)
ax.set_title("CLs vs Signal Strength", fontsize=16)
ax.set_xlim(0, min(obs_limit * 3, poi_values[-1]))
ax.set_ylim(0, 1)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("limit_cls_scan.png", dpi=150)
print("  Saved: limit_cls_scan.png")
plt.show()

print("\nDone!")
