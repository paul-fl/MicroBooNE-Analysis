#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# Paths
base_dir = os.path.expanduser("~/Msci/outputs/inference")
steps = ["7346", "8441", "9821"]
signal_file = "dm_signal_only_test_set_DM-CNN_scores_inference_DM_CNN_pf522.csv"
background_file = "cosmics_corsika_test_set_DM-CNN_scores_inference_DM_CNN_pf522.csv"
plot_outdir = os.path.expanduser("~/Msci/plots/inference/ROC_curves")
os.makedirs(plot_outdir, exist_ok=True)

# Style (MicroBooNE-like)
plt.rcParams.update({
    "font.size": 11,
    "font.family": "sans-serif",
    "axes.linewidth": 1.2,
    "axes.labelpad": 8,
    "axes.titlesize": 13,
    "legend.framealpha": 1,
    "legend.edgecolor": "none",
})

fig, ax = plt.subplots(figsize=(4.6, 4.6))

for step in steps:
    step_dir = os.path.join(base_dir, f"step_{step}")
    sig_path = os.path.join(step_dir, signal_file)
    bkg_path = os.path.join(step_dir, background_file)
    if not (os.path.exists(sig_path) and os.path.exists(bkg_path)):
        continue

    df_sig = pd.read_csv(sig_path)
    df_bkg = pd.read_csv(bkg_path)
    df_sig["label"] = 1
    df_bkg["label"] = 0
    df = pd.concat([df_sig, df_bkg], ignore_index=True)

    fpr, tpr, _ = roc_curve(df["label"], df["signal_score"])
    roc_auc = auc(tpr, 1 - fpr)
    ax.plot(
        tpr,
        1 - fpr,
        lw=1.8,                      # thinner lines
        label=f"Step {step} [AUC = {roc_auc:.3f}]",
    )

# Axes settings
ax.set_xlabel("Signal efficiency")
ax.set_ylabel("Background rejection efficiency")

# Extend both axes beyond 1.0
ax.set_xlim(0.0, 1.05)
ax.set_ylim(0.0, 1.05)

# Add clean ticks
ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

# Clean legend and title
ax.legend(loc="lower left", fontsize=9, frameon=True, facecolor="white")
ax.set_title("DM-CNN ROC Curves", pad=8)
ax.set_aspect("equal", adjustable="box")

fig.tight_layout()

output_path = os.path.join(plot_outdir, "DM_CNN_ROC_refined.png")
fig.savefig(output_path, dpi=300)
print(f"✅ Saved to {output_path}")
