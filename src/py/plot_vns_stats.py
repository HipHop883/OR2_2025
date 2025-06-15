import pandas as pd
import matplotlib.pyplot as plt

# Load VNS log data
df = pd.read_csv("../../logs/vns_log.csv")

# Identify when best_cost is improved
best_updates = df[df['best_cost'].diff() < 0]

# Setup figure
fig, ax1 = plt.subplots(figsize=(14, 6))

# Plot current and best cost
ax1.plot(df['iteration'], df['current_cost'], label='Current Cost', color='blue', linewidth=1)
ax1.plot(df['iteration'], df['best_cost'], label='Best Cost (Incumbent)', linestyle='--', color='orange', linewidth=1.5)

# Highlight new bests
ax1.scatter(best_updates['iteration'], best_updates['best_cost'],
            color='red', marker='x', s=40, label='New Best Found')

# Configure primary axis
ax1.set_xlabel('Iteration', fontsize=12)
ax1.set_ylabel('Solution Cost', fontsize=12, color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True, linestyle='--', alpha=0.5)

# Add secondary axis for 3-opt kicks
ax2 = ax1.twinx()
ax2.plot(df['iteration'], df['kicks'], label='3-opt Kicks',
         color='green', linestyle='-.', linewidth=1.3, alpha=0.6)
ax2.set_ylabel('3-opt Kicks', fontsize=12, color='green')
ax2.tick_params(axis='y', labelcolor='green')

# Title and legend
plt.title("VNS on att48: Solution Cost and 3-opt Kicks vs Iteration", fontsize=14)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)

# Layout and display
plt.tight_layout()
plt.show()
