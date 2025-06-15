import pandas as pd
import matplotlib.pyplot as plt

# Load log data
df = pd.read_csv("../../logs/tabu_log.csv")

# Detect best updates
best_updates = df[df['best_cost'].diff() < 0]

# Detect diversification spikes (optional, based on cost spikes)
# This is a heuristic: large positive jumps in cost
div_threshold = 500
diversifications = df[df['current_cost'].diff() > div_threshold]

# Setup figure
fig, ax1 = plt.subplots(figsize=(14, 6))

# Plot costs
ax1.plot(df['iteration'], df['current_cost'], label='Current Cost', color='blue', linewidth=1)
ax1.plot(df['iteration'], df['best_cost'], label='Best Cost (Incumbent)', linestyle='--', color='orange', linewidth=1.5)

# Plot new bests
ax1.scatter(best_updates['iteration'], best_updates['best_cost'],
            color='red', marker='x', s=40, label='New Best Found')

# Optional: draw vertical lines for diversification (spikes in cost)
for _, row in diversifications.iterrows():
    ax1.axvline(row['iteration'], color='gray', linestyle=':', alpha=0.3)

# Left axis (Cost)
ax1.set_xlabel('Iteration', fontsize=12)
ax1.set_ylabel('Cost', fontsize=12, color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True, linestyle='--', alpha=0.5)

# Right axis (Tenure)
ax2 = ax1.twinx()
ax2.plot(df['iteration'], df['tenure'], label='Tabu Tenure',
         color='purple', linestyle='-.', alpha=0.6, linewidth=1.5)
ax2.set_ylabel('Tabu Tenure', fontsize=12, color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# Title and legend
plt.title("Tabu Search on att48: Cost, Best Found, and Tenure vs Iteration", fontsize=14)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)

# Layout and export
plt.tight_layout()
plt.show()
