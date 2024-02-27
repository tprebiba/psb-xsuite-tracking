#%%
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
from simulation_parameters import parameters as p


#%%
output = np.load('../output/emittances.npy')
#output = np.load('/home/tprebiba/afs/temp/emittances.npy')
macrosize = p['bunch_intensity']/p['n_part']
intensity = output[:,0]*macrosize
turns = np.arange(len(intensity))
epsx = output[:,1]
epsy = output[:,2]
epsz = output[:,3]


# %%
f, axs = plt.subplots(2,1,figsize=(6,6), sharex=True)
fontsize=15
ax = axs[0]
ax.tick_params(axis='both', which='major', labelsize=fontsize)
ax.set_ylabel('Intensity [10$^{10}$ ppb]', fontsize=fontsize)
#ax.set_xlabel('Turn', fontsize=fontsize)
ax.plot(turns, intensity/1e10, color='blue')
ax = axs[1]
ax.tick_params(axis='both', which='major', labelsize=fontsize)
ax.set_ylabel('Emittance [mm mrad]', fontsize=fontsize)
ax.set_xlabel('Turn', fontsize=fontsize)
ax.set_ylim(0.2, 0.38)
ax.plot(turns, epsx*1e6, color='green', label='x')
ax.plot(turns, epsy*1e6, color='red', label='y')
#ax.plot(turns, epsz, color='green', label='z')
ax.legend(loc=0, fontsize=fontsize)
#ax.set_xlim(0,1000)
plt.tight_layout()
plt.show()
# %%
