import numpy as np
import pandas as pd
import pyhf
import matplotlib.pyplot as plt

# -----------------------------
# USER SETTINGS
# -----------------------------
nominal_eps = 1e-3
alpha = 0.1
ratio = "0.6"
tag = "CNN"

fraction_outside = 1.1
flux_uncert = 0.22
pot_uncert = 0.02
flat_uncertainty = False
uncertainty = 0.5

tests = False
masses = ["0.05"]   # Just like Luis had for single-run checks

# ===================================================================
# STEP 1 — LOAD *YOUR* CSV SCORES AND BUILD HISTOGRAMS
# ===================================================================

# File paths
RUN1 = "/home/paul/Msci/outputs/inference/run1/"
RUN3 = "/home/paul/Msci/outputs/inference/run3/"

def load_scores(path, fname):
    df = pd.read_csv(path + fname)
    return df["signal_score"].values

# --- Run 1 CSVs ---
sig1  = load_scores(RUN1, "run1_dt_ratio_0.6_ma_0.05_pi0_larcv_cropped_scores.csv")
nu1   = load_scores(RUN1, "run1_NuMI_nu_overlay_larcv_cropped_scores.csv")
dirt1 = load_scores(RUN1, "run1_NuMI_dirt_larcv_cropped_scores.csv")
off1  = load_scores(RUN1, "run1_offbeam_larcv_cropped_full_set_scores.csv")

# --- Run 3 CSVs ---
sig3  = load_scores(RUN3, "run3_dt_ratio_0.6_ma_0.05_pi0_larcv_cropped_scores.csv")
nu3   = load_scores(RUN3, "run3_nu_overlay_larcv_cropped_scores.csv")
dirt3 = load_scores(RUN3, "run3_dirt_larcv_cropped_scores.csv")
off3  = load_scores(RUN3, "run3_offbeam_larcv_cropped_scores.csv")

# Luis uses ~20 bins (0–1)
bins = np.linspace(0, 1, 21)

# Table 7.4 weights
w1 = {"nu":0.09, "dirt":0.10, "off":0.66}
w3 = {"nu":0.24, "dirt":0.17, "off":0.30}

def weighted_hist(scores, w):
    return np.histogram(scores, bins=bins, weights=np.ones_like(scores)*w)[0]

# ----- RUN 1 histograms -----
sig1_h  = np.histogram(sig1, bins=bins)[0]
nu1_h   = weighted_hist(nu1,   w1["nu"])
dirt1_h = weighted_hist(dirt1, w1["dirt"])
off1_h  = weighted_hist(off1,  w1["off"])
bkg1_h  = nu1_h + dirt1_h + off1_h
data1_h = bkg1_h  # no real data → Asimov

# ----- RUN 3 histograms -----
sig3_h  = np.histogram(sig3, bins=bins)[0]
nu3_h   = weighted_hist(nu3,   w3["nu"])
dirt3_h = weighted_hist(dirt3, w3["dirt"])
off3_h  = weighted_hist(off3,  w3["off"])
bkg3_h  = nu3_h + dirt3_h + off3_h
data3_h = bkg3_h

# ----- COMBINE RUNS -----
signal_hist = np.concatenate([sig1_h, sig3_h])
background_hist = np.concatenate([bkg1_h, bkg3_h])
data_hist = np.concatenate([data1_h, data3_h])

# Luis’s background uncertainties (shapesys)
back_sigma = np.sqrt(background_hist)  # minimal approximation

print("Loaded and histogrammed CNN scores.")
print("Total signal events =", signal_hist.sum())
print("Total background events =", background_hist.sum())
print("S/B =", signal_hist.sum()/background_hist.sum())


# ===================================================================
# STEP 2 — LUIS'S MODEL (UNCHANGED)
# ===================================================================

for mass in masses:

    # SIGNAL SHAPE UNCERTAINTY: use sqrt(signal) as approx
    sigma_sig = np.sqrt(signal_hist)

    model = pyhf.Model({
        "channels": [{
            "name": "singlechannel",
            "samples": [
                {
                    "name": "signal",
                    "data": signal_hist.tolist(),
                    "modifiers": [
                        {"name":"mu", "type":"normfactor", "data":None},
                        {"name":"uncorr_siguncrt", "type":"shapesys", "data": sigma_sig.tolist()},
                        {"name":"pot_correlated", "type":"normsys", "data":{"hi":1+pot_uncert, "lo":1-pot_uncert}},
                        {"name":"flux_correlated","type":"normsys","data":{"hi":1+flux_uncert, "lo":1-flux_uncert}},
                    ]
                },
                {
                    "name": "background",
                    "data": background_hist.tolist(),
                    "modifiers": [
                        {"name":"uncorr_bkguncrt", "type":"shapesys", "data": back_sigma.tolist()},
                        {"name":"pot_correlated", "type":"normsys","data":{"hi":1+pot_uncert, "lo":1-pot_uncert}},
                    ]
                }
            ]
        }]
    })

    obs = data_hist.tolist() + model.config.auxdata

    poi_values = np.linspace(0., 10., 100)
    obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upperlimit(
        obs, model, poi_values, level=0.1, return_results=True
    )

    print("\n===== RESULTS =====")
    print(f"Observed μ95 = {obs_limit:.4f}")
    print(f"Expected μ95 = {exp_limits[2]:.4f}")

    epsilon2_obs = (nominal_eps**2)*np.sqrt(obs_limit)
    print(f"Observed ε² = {epsilon2_obs:.3e}")


# ===================================================================
# STEP 3 — LUIS'S BRAZIL PLOT (UNCHANGED)
# ===================================================================

from pyhf.contrib.viz import brazil

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)

brazil.plot_results(poi_values, results, test_size=0.1, ax=ax)

ax.set_xlabel("μ")
ax.set_ylabel("CLs")
ax.set_title("CLs — Luis-Style Limit (CNN Histograms)")

plt.tight_layout()
plt.show()

