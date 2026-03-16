#!/usr/bin/env python3
"""Plot MPI scaling curves (speedup & efficiency) for python-svof RT benchmarks."""

import matplotlib.pyplot as plt
import numpy as np

# --- Data from scaling runs (3/15/26) ---

data = {
    'RT 160x200 (32k cells, t=6.0)': {
        'np':   [1,      2,      4,      8,      12,     16,     20,     24],
        'time': [97.1554, 52.5382, 28.8253, 21.5877, 19.7873, 18.9986, 18.5512, 21.4087],
    },
    'RT 320x400 (128k cells, t=3.0)': {
        'np':   [1,       2,       4,       8,       16,      18,      20,      24],
        'time': [476.8866, 177.2468, 99.1557, 66.8103, 41.9209, 36.9211, 34.9755, 40.4455],
    },
    'RT 640x800 (512k cells, t=1.0)': {
        'np':   [1,        2,       4,       8,       16,      18,      20,      24],
        'time': [1768.575, 525.154, 299.635, 168.5934, 86.4175, 77.6180, 72.3584, 76.3010],
    },
}

markers = ['o', 's', '^']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

for i, (label, d) in enumerate(data.items()):
    nps = np.array(d['np'])
    times = np.array(d['time'])
    t1 = times[0]
    speedup = t1 / times
    efficiency = speedup / nps * 100

    ax1.plot(nps, speedup, marker=markers[i], color=colors[i],
             linewidth=2, markersize=7, label=label)
    ax2.plot(nps, efficiency, marker=markers[i], color=colors[i],
             linewidth=2, markersize=7, label=label)

# Ideal scaling line
np_max = 24
ax1.plot([1, np_max], [1, np_max], 'k--', linewidth=1, alpha=0.5, label='Ideal')

ax1.set_xlabel('MPI ranks', fontsize=12)
ax1.set_ylabel('Speedup (T₁ / Tₙ)', fontsize=12)
ax1.set_title('Strong Scaling — Speedup', fontsize=13)
ax1.legend(fontsize=9, loc='upper left')
ax1.set_xlim(0, np_max + 1)
ax1.set_ylim(0, np_max + 1)
ax1.set_xticks([1, 2, 4, 8, 12, 16, 20, 24])
ax1.grid(True, alpha=0.3)

ax2.axhline(100, color='k', linestyle='--', linewidth=1, alpha=0.5, label='Ideal')
ax2.set_xlabel('MPI ranks', fontsize=12)
ax2.set_ylabel('Parallel efficiency (%)', fontsize=12)
ax2.set_title('Strong Scaling — Efficiency', fontsize=13)
ax2.legend(fontsize=9, loc='upper right')
ax2.set_xlim(0, np_max + 1)
ax2.set_ylim(0, 175)
ax2.set_xticks([1, 2, 4, 8, 12, 16, 20, 24])
ax2.grid(True, alpha=0.3)

fig.suptitle('python-svof MPI Strong Scaling — Rayleigh-Taylor', fontsize=14,
             fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig('scaling_curves.png', dpi=150, bbox_inches='tight')
fig.savefig('scaling_curves.pdf', bbox_inches='tight')
print("Saved scaling_curves.png and scaling_curves.pdf")
plt.show()
