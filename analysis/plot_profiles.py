#%%
import numpy as np
import matplotlib.pyplot as plt
import json
import pandas as pd
import matplotlib.cm as cm
import beamtools.fitting as btf


# %%
###########################################
# Load Xsuite data
###########################################
source_dir = './'
directories = [
    'lhc10',
    'lhc15',
    'lhc20',
    'lhc25',
    'lhc30',
    'lhc35'
]
prof = {}
for dir in directories:
    #file = source_dir + dir+'/output/particles_turn_24999.json'
    file = source_dir + dir+'/input/particles_initial.json'
    df = pd.read_json(file)
    y = df['y']
    yprof, ypos = np.histogram(y*1e3, bins=200, range=(-29,29))
    prof[dir] = {'yprof': yprof, 'ypos': ypos[0:-1]}
colors = cm.rainbow(np.linspace(0, 1, len(directories)))


# %%
###########################################
# Plot intensity and emittance vs turn
###########################################
# profile
f,ax = plt.subplots(1,1,figsize=(6,4.5),facecolor='white')
font_size = 22
ax.set_xlabel('y [mm]', fontsize=font_size)
ax.tick_params(axis='y',which='both',left=False,labelleft=False)
ax.tick_params(labelsize=font_size)

dir = 'lhc10'
x = prof[dir]['ypos']
y = prof[dir]['yprof']
popt,pcov,popte,x_gaussianfit,y_gaussianfit = btf.makeGaussianFit(x,y)
q0 = 1.1
p0 = [popt[1], q0,  1/popt[2]**2/(5-3*q0), popt[0]] 
poptq, pcovq, poptqe, x_qgaussianfit, y_qgaussianfit = btf.makeQGaussianFit(x,y,p0)
print('q = ' ,poptq[1])
ax.plot(x,y, '.', lw=1., color='black', label=r'$N_b \approx %s0$, $q=%1.2f$'%(dir.split('lhc')[1], poptq[1]))
ax.plot(x_gaussianfit, y_gaussianfit, color='green')
ax.plot(x_qgaussianfit, y_qgaussianfit, color='red')

ax.set_xlim(-11,11)
#ax.set_yscale('log')
#ax.set_xlim(-27,27)
#ax.set_ylim(100)
ax.legend(loc=0, fontsize=font_size-3, numpoints=1)
f.tight_layout()


# %%
