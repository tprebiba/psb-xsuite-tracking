#%%
import numpy as np
import pylab as plt
from scipy.interpolate import griddata
from PyNAFF import naff
from resonance_lines import resonance_lines
import json


#%%
x=np.load('../output/x.npy')
y=np.load('../output/y.npy')
p = {}
with open('../tables/tunes.json', 'r') as fid:
     d = json.load(fid)
for key in d:
     p[key] = d[key]
try:
    qx_pic = np.load('../output/qx_pic.npy')
    qy_pic = np.load('../output/qy_pic.npy')
except:
    qx_pic = np.nan
    qy_pic = np.nan
    print('No pic data')


#%%
Qx1, Qx = [], []
Qy1, Qy = [], []

for i in range(len(x)):
    try:
        Qx1.append(naff(x[i]-np.mean(x[i]), turns=600)[0][1])
    except:
        Qx1.append(np.nan)
    try:
        Qy1.append(naff(y[i]-np.mean(y[i]), turns=600)[0][1])
    except:
        Qy1.append(np.nan)
    try:
        Qx.append(naff(x[i]-np.mean(x[i]), skipTurns=599, turns=600)[0][1])
    except:
        Qx.append(np.nan)
    try:
        Qy.append(naff(y[i]-np.mean(y[i]), skipTurns=599, turns=600)[0][1])
    except:
        Qy.append(np.nan)
Qx1=4.0+np.array(Qx1)
Qx=4.0+np.array(Qx)
Qy1=4.0+np.array(Qy1)
Qy=4.0+np.array(Qy)+0.1
d = np.log(np.sqrt( (Qx-Qx1)**2 + (Qy-Qy1)**2 ))


#%%
r=resonance_lines([4.0,4.3],[4.4,4.7], [1,2,3,4], 16)
#r=resonance_lines([4.06,4.3],[4.06,4.3], [1,2,3,4], 16)
fontsize=17
f, ax = plt.subplots(1,figsize=(6,6))
my_cmap = plt.cm.jet
my_cmap.set_under('w',1)

qxnew = np.copy(qx_pic)#[mask]
qynew = np.copy(qy_pic)#[mask]
qxnew[np.where(qxnew<0)] += 1
qynew[np.where(qynew<0)] += 1
ax.hist2d(qxnew+4, qynew+4,
          bins=500,#400, 
          cmap = my_cmap, 
          vmin=10, #5, 
          range=[[r.Qx_min, r.Qx_max], [r.Qy_min, r.Qy_max]])#, norm=mcolors.PowerNorm(gamma))
#ax.plot(qx_pic+4, qy_pic+4,'k.',ms=0.5)
#ax.plot(qxnew+4, qynew+4,'k.',ms=0.5)

ax.scatter(Qx, Qy,4, d, 'o',lw = 0.1,zorder=10, cmap=my_cmap)
ax.plot([p['qx_ini']],[p['qy_ini']],'ko',zorder=1e5)
ax.set_xlabel('$\mathrm{Q_x}$', fontsize=fontsize)
ax.set_ylabel('$\mathrm{Q_y}$', fontsize=fontsize)
ax.tick_params(axis='both', labelsize=fontsize)
#plt.clim(-20.5,-4.5)

r.plot_resonance(f)

ax.xaxis.label.set_size(fontsize)
ax.yaxis.label.set_size(fontsize)
ax.tick_params(labelsize=fontsize)
#cbar=plt.colorbar()
#cbar.set_label('d',fontsize='18')
#cbar.ax.tick_params(labelsize='18')
plt.tight_layout()


#%%
fig2=plt.figure()
XX,YY = np.meshgrid(np.unique(x[:,0]), np.unique(y[:,0]))
Z = griddata((x[:,0],y[:,0]), d, (XX,YY), method='linear')
Zm = np.ma.masked_invalid(Z)
fig2.suptitle('Initial Distribution', fontsize='20')
plt.pcolormesh(XX,YY,Zm,cmap=plt.cm.jet)
# plt.scatter(x[:,0],y[:,0],4, d, 'o',lw = 0.1,zorder=10, cmap=plt.cm.jet)
plt.tick_params(axis='both', labelsize='18')
plt.xlabel('x [m]', fontsize='20')
plt.ylabel('y [m]', fontsize='20')
plt.clim(-20.5,-4.5)
cbar=plt.colorbar()
cbar.set_label('d',fontsize='18')
cbar.ax.tick_params(labelsize='18')
#fig2.savefig('test_initial_distribution.png')
plt.show()
# %%
