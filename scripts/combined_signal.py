import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------
# Paths
# -----------------------------
data_dir = "/home/paul/Msci/outputs/inference/run3"
plot_dir = "/home/paul/Msci/plots/inference/run3"
os.makedirs(plot_dir, exist_ok=True)

# -----------------------------
# Background definitions (single file each)
# -----------------------------
# filename, label, color, weight
backgrounds = [
    ("run3_nu_overlay_larcv_cropped_scores.csv", "In cryo", "#8f0ad1", 0.24),
    ("run3_dirt_larcv_cropped_scores.csv", "Out of cryo", "#1d0265", 0.17),
    ("run3_offbeam_larcv_cropped_scores.csv", "Beam-off", "#35b6e3", 0.30),
]

signal = ("run3_dt_ratio_0.6_ma_0.05_pi0_larcv_cropped_scores.csv", "Dark Trident", "red")

bins = np.linspace(0, 1, 15)

bg_scores = []
bg_weights = []
bg_labels = []
bg_colors = []

# -----------------------------
# Load background samples
# -----------------------------
for filename, label, color, weight in backgrounds:
    filepath = os.path.join(data_dir, filename)
    if not os.path.exists(filepath):
        print(f"Warning: File not found: {filepath}")
        continue
    
    df = pd.read_csv(filepath)
    scores = df["signal_score"].to_numpy()
    scores = scores[(scores >= 0) & (scores <= 1)]
    
    bg_scores.append(scores)
    bg_weights.append(np.ones_like(scores) * weight)
    bg_labels.append(label)
    bg_colors.append(color)
    
    print(f"Loaded {label}: {len(scores)} events")

# -----------------------------
# Load signal
# -----------------------------
sig_filepath = os.path.join(data_dir, signal[0])
if not os.path.exists(sig_filepath):
    raise RuntimeError(f"Signal file not found: {sig_filepath}")

df_sig = pd.read_csv(sig_filepath)
sig_scores = df_sig["signal_score"].to_numpy()
sig_scores = sig_scores[(sig_scores >= 0) & (sig_scores <= 1)]

print(f"Loaded {signal[1]}: {len(sig_scores)} events")

# -----------------------------
# PLOT 1: Standard CNN score histogram (0–1)
# -----------------------------
plt.figure(figsize=(8, 6))

plt.hist(
    bg_scores,
    bins=bins,
    stacked=True,
    weights=bg_weights,
    color=bg_colors,
    label=bg_labels,
    edgecolor="black",
    linewidth=0.3,
)

plt.hist(
    sig_scores,
    bins=bins,
    histtype="step",
    color=signal[2],
    linewidth=2,
    label=signal[1],
)

plt.axvline(0.5, color="black", linestyle="--", linewidth=1)
plt.xlim(0, 1)
plt.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.7)
plt.xlabel("CNN Signal Score", fontsize=12)
plt.ylabel("Events", fontsize=12)
plt.legend(framealpha=1, edgecolor="black", fontsize=10)
plt.tight_layout()

out1 = os.path.join(plot_dir, "run1_signal_background_stacked.png")
plt.savefig(out1, dpi=300)
print(f"Saved: {out1}")

# -----------------------------
# LOGIT PLOT
# -----------------------------
def logit(x):
    eps = 1e-12
    x = np.clip(x, eps, 1 - eps)
    return np.log(x / (1 - x))

bg_sr = [s[s > 0.5] for s in bg_scores]
bg_sr_weights = [w[s > 0.5] for w, s in zip(bg_weights, bg_scores)]
sig_sr = sig_scores[sig_scores > 0.5]

bg_logit = [logit(s) for s in bg_sr]
sig_logit = logit(sig_sr)

logit_bins = np.array([0, 1.4, 2.8, 4.1, 5.5, 6.9, 8])
signal_scale = 0.075

plt.figure(figsize=(8, 6))

plt.hist(
    bg_logit,
    bins=logit_bins,
    stacked=True,
    weights=bg_sr_weights,
    color=bg_colors,
    edgecolor="black",
    linewidth=0.3,
    label=bg_labels,
)

plt.hist(
    sig_logit,
    bins=logit_bins,
    histtype="step",
    weights=np.ones_like(sig_logit) * signal_scale,
    color="red",
    linewidth=2,
    label="Dark Trident (scaled)",
)

plt.xlabel("CNN Signal Score (logit)")
plt.ylabel("Events")
plt.xlim(0, 7)
plt.legend(framealpha=1, edgecolor="black")
plt.tight_layout()

out2 = os.path.join(plot_dir, "run1_signal_background_logit.png")
plt.savefig(out2, dpi=300)
print(f"Saved logit: {out2}")

plt.show()
