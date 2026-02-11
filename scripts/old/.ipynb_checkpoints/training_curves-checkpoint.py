import pandas as pd
import matplotlib.pyplot as plt

# Load the training metrics
df = pd.read_csv("/home/paul/Msci/metrics/DM-CNN_training_metrics_20251208-03_21_PM_multiclass_run.csv")

# Create figure with 2 subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Loss curves
axes[0].plot(df['step'], df['train_loss'], label='Train Loss', alpha=0.7)
axes[0].plot(df['step'], df['test_loss'], label='Test Loss', alpha=0.7)
axes[0].set_xlabel('Step')
axes[0].set_ylabel('Loss')
axes[0].set_title('Training and Test Loss')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Plot 2: Accuracy curves
axes[1].plot(df['step'], df['train_accu'], label='Train Accuracy', alpha=0.7)
axes[1].plot(df['step'], df['test_accu'], label='Test Accuracy', alpha=0.7)
axes[1].set_xlabel('Step')
axes[1].set_ylabel('Accuracy')
axes[1].set_title('Training and Test Accuracy')
axes[1].legend()
axes[1].grid(True, alpha=0.3)
axes[1].set_ylim(0, 1.05)

plt.tight_layout()
plt.savefig("/home/paul/Msci/plots/multiclass_training_curves.png", dpi=300)
plt.savefig("/home/paul/Msci/plots/multiclass_training_curves.pdf")
plt.show()

# Print final stats
print("\n=== Final Training Stats ===")
print(f"Final train accuracy: {df['train_accu'].iloc[-1]:.3f}")
print(f"Final test accuracy: {df['test_accu'].iloc[-1]:.3f}")
print(f"Final train loss: {df['train_loss'].iloc[-1]:.3f}")
print(f"Final test loss: {df['test_loss'].iloc[-1]:.3f}")
print(f"Total steps: {df['step'].iloc[-1]:.0f}")
print(f"Total epochs: {df['epoch'].iloc[-1]:.0f}")

# Check if still improving
last_100 = df.tail(10)
print(f"\nLast 10 steps - avg test loss: {last_100['test_loss'].mean():.3f}")
first_100 = df.iloc[50:60]  # Skip initial instability
print(f"Steps 500-600 - avg test loss: {first_100['test_loss'].mean():.3f}")
