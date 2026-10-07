import matplotlib.pyplot as plt
import numpy as np
 
cycles = [0, 1, 2, 3, 4]
num_WTs = [53700, 3210, 158, 46, 5]
num_pBADs = [3670000, 1370, 31, 2, 1]
 
# Survival fraction relative to cycle 0
survival_WT = [n / num_WTs[0] for n in num_WTs]
survival_pBAD = [n / num_pBADs[0] for n in num_pBADs]
 
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
 
# Panel 1: raw colony counts, log y-axis
ax = axes[0]
ax.plot(cycles, num_WTs, 'o-', label='WT', color='#2A6F97')
ax.plot(cycles, num_pBADs, 's-', label='pBAD (empty vector)', color='#B3541E')
ax.set_yscale('log')
ax.set_xlabel('Freeze-thaw cycle')
ax.set_ylabel('Colony count (CFU)')
ax.set_title('Raw colony counts')
ax.legend()
ax.grid(True, which='both', alpha=0.3)
 
# Panel 2: log10 survival fraction relative to cycle 0
ax = axes[1]
ax.plot(cycles, np.log10(survival_WT), 'o-', label='WT', color='#2A6F97')
ax.plot(cycles, np.log10(survival_pBAD), 's-', label='pBAD (empty vector)', color='#B3541E')
ax.set_xlabel('Freeze-thaw cycle')
ax.set_ylabel(r'$\log_{10}(N/N_0)$')
ax.set_title('Log survival relative to cycle 0')
ax.legend()
ax.grid(True, alpha=0.3)
 
plt.tight_layout()
plt.savefig('freeze_thaw_survival.png', dpi=300, bbox_inches='tight')
plt.show()
