#%%
import numpy as np
import pylab as plt
from scipy.interpolate import griddata
from PyNAFF import naff
from resonance_lines import resonance_lines
import matplotlib.gridspec as gridspec
import json

#%%
qx_pic = np.load('../output/qx_pic.npy')
qy_pic = np.load('../output/qy_pic.npy')
qxnew = np.copy(qx_pic)
qynew = np.copy(qy_pic)
qxnew[np.where(qxnew<0)] += 1
qynew[np.where(qynew<0)] += 1
qxnew = qxnew+4
qynew = qynew+4

#%%
my_cmap = plt.cm.jet
my_cmap.set_under('w',1)
fontsize=20
bins = 1000
vmin = 10
#vmin = 1
r=resonance_lines([4.1,4.2],[4.5,4.7], [1,2,3,4], 16)

fig = plt.figure(figsize=(10,8))
gs = gridspec.GridSpec(6, 6)
ax = plt.subplot(gs[1:6, 0:4])
ax_xDist = plt.subplot(gs[0, 0:4],sharex=ax)
ax_yDist = plt.subplot(gs[1:6, 4],sharey=ax)

ax.hist2d(qxnew, qynew,
          bins=500, cmap = my_cmap, vmin=vmin, 
          range=[[r.Qx_min, r.Qx_max], [r.Qy_min, r.Qy_max]]#, norm=mcolors.PowerNorm(gamma)
          )
ax.plot(4.15, 4.65, '*', ms=20, color='k', zorder=1e5)

#r.plot_resonance(ax)
ax.set_xlabel('$\mathrm{Q_x}$', fontsize=fontsize)
ax.set_ylabel('$\mathrm{Q_y}$', fontsize=fontsize)
ax.tick_params(axis='both', labelsize=fontsize)
ax.set_xlim(4.0, 4.2)
ax.set_ylim(4.5, 4.7)
ax.axhline(4.45, color='k', ls='--')

ax_xDist.hist(qxnew,bins=bins,align='mid', color='black')
ax_xDist.tick_params(axis='x', labelleft=False, labelright=False, labeltop=False, labelbottom=False)
ax_xDist.tick_params(axis='y', labelleft=False, labelright=False, labeltop=False, labelbottom=False)
nsigmas = 2.5
ax_xDist.axvline(np.mean(qxnew)-nsigmas*np.std(qxnew), color='red', ls='-')
ax_xDist.axvline(np.mean(qxnew)+nsigmas*np.std(qxnew), color='red', ls='-')
print('Total tune spread in x: ', nsigmas*np.std(qxnew)*2)
#ax_xDist.set(ylabel='count')
#ax_xCumDist = ax_xDist.twinx()
#ax_xCumDist.hist(qxnew,bins=bins,cumulative=True,histtype='step',density=True,color='r',align='mid')
#ax_xCumDist.tick_params('y', colors='r')
#ax_xCumDist.set_ylabel('cumulative',color='r')

ax_yDist.hist(qynew,bins=bins,orientation='horizontal',align='mid', color='black')
ax_yDist.tick_params(axis='x', labelleft=False, labelright=False, labeltop=False, labelbottom=False)
ax_yDist.tick_params(axis='y', labelleft=False, labelright=False, labeltop=False, labelbottom=False)
nsigmas = 2.5
ax_yDist.axhline(np.mean(qynew)-nsigmas*np.std(qynew), color='red', ls='-')
ax_yDist.axhline(np.mean(qynew)+nsigmas*np.std(qynew), color='red', ls='-')
print('Total tune spread in y: ', nsigmas*np.std(qynew)*2)
#ax_yDist.set(xlabel='count')
#ax_yCumDist = ax_yDist.twiny()
#ax_yCumDist.hist(qynew,bins=bins,cumulative=True,histtype='step',density=True,color='r',align='mid',orientation='horizontal')
#ax_yCumDist.tick_params('x', colors='r')
#ax_yCumDist.set_xlabel('cumulative',color='r')

fig.tight_layout()
# %%
