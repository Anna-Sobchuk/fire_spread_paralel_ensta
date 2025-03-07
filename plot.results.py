import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('results.csv')
baseline_total = df[df.Threads == 1].TotalTime.values[0]
baseline_update = df[df.Threads == 1].UpdateTime.values[0]
df['SpeedupTotal'] = baseline_total / df.TotalTime
df['SpeedupUpdate'] = baseline_update / df.UpdateTime

plt.figure(figsize=(8, 5))
x = df.Threads
xticks = df.Threads.unique()


plt.plot(x, df.SpeedupTotal, 'o--', label='Total Acceleration', markersize=8)
plt.plot(x, df.SpeedupUpdate, 's-', label='TotUpdateStepal Acceleration', markersize=8)

# Formatting
plt.title(f"Parallel Acceleration")
plt.xlabel('Number of Threads')
plt.ylabel('Speedup Factor')
plt.xticks(xticks)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()

# Save and show
plt.savefig('acceleration_plot.png', dpi=300)
plt.show()