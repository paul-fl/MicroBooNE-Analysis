import pandas as pd

# Load inference results
df = pd.read_csv("../plots/inference/MPID_test_set_full_DM-CNN_scores_inference_DM_CNN_pf522.csv")

# Define categories
uncertain = df[(df["signal_score"] > 0.4) & (df["signal_score"] < 0.6)]
confident_signal = df[df["signal_score"] > 0.9]
confident_background = df[df["signal_score"] < 0.1]

print("\n Uncertain entries (signal_score ~0.5):")
print(uncertain[["entry_number", "signal_score"]].head(10).to_string(index=False))

print("\n Confident signal entries (signal_score > 0.9):")
print(confident_signal[["entry_number", "signal_score"]].head(10).to_string(index=False))

print("\n Confident background entries (signal_score < 0.1):")
print(confident_background[["entry_number", "signal_score"]].head(10).to_string(index=False))

