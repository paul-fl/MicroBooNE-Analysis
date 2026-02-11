import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Path where your CSVs are stored
DATA_DIR = "/home/paul/Msci/dm_sets_multiclass/run1_signal"
OUTPUT_DIR = "/home/paul/Msci/plots/inference_multiclass/run1"

# Map files to their labels (update paths to match your folder structure)
file_mapping = [
    ("run1_NuMI_nu_overlay_larcv_cropped_scores_multiclass.csv", "In cryo"),
    ("run1_NuMI_dirt_larcv_cropped_scores_multiclass.csv", "Out of cryo"),
    ("run1_offbeam_larcv_cropped_full_set_scores_multiclass.csv", "Beam-off"),
    ("run1_dt_ratio_0.6_ma_0.05_pi0_larcv_cropped_scores_multiclass.csv", "Dark Trident"),
]
CLASS_NAMES = ["NCπ0", "Cosmics", "DM Signal"]

# Load data
dfs = []
for filename, label in file_mapping:
    path = os.path.join(DATA_DIR, filename)
    df = pd.read_csv(path)
    df["true_label"] = label
    dfs.append(df)

full_df = pd.concat(dfs, ignore_index=True)

# Get predicted class
if 'pred_class' not in full_df.columns:
    prob_cols = ['p0', 'p1', 'p2']
    full_df['pred_class'] = full_df[prob_cols].values.argmax(axis=1)

# Build confusion matrix
true_labels = ["In cryo", "Out of cryo", "Beam-off", "Dark Trident"]
cm = np.zeros((4, 3))

for i, true_label in enumerate(true_labels):
    subset = full_df[full_df['true_label'] == true_label]
    for j in range(3):
        cm[i, j] = (subset['pred_class'] == j).sum()

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(cm, cmap="Blues")
plt.colorbar(im, ax=ax)

ax.set_xticks(range(3))
ax.set_yticks(range(4))
ax.set_xticklabels(CLASS_NAMES, rotation=45, ha="right")
ax.set_yticklabels(true_labels)
ax.set_xlabel("Predicted Class")
ax.set_ylabel("Actual File Type")
ax.set_title("Run1 Multiclass Confusion Matrix (3 classes)")

# Annotate
for i in range(4):
    for j in range(3):
        ax.text(j, i, int(cm[i, j]), ha="center", va="center")

plt.tight_layout()
os.makedirs(OUTPUT_DIR, exist_ok=True)
plt.savefig(os.path.join(OUTPUT_DIR, "run1_confusion_matrix_3class.png"), dpi=300)
plt.savefig(os.path.join(OUTPUT_DIR, "run1_confusion_matrix_3class.pdf"))
plt.show()
