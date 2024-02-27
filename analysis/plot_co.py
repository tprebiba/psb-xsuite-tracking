import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable


f, ax = plt.subplots(1,1,figsize=(9,6))
fontsize=20
ax.set_xlabel('s [m]', fontsize=fontsize)
ax.set_ylabel('x [mm]', fontsize=fontsize)
ax.tick_params(labelsize=fontsize)

toplot = range(0,1000,9)
colors = cm.rainbow(np.linspace(0,1,len(toplot)))
for c,i in zip(colors,toplot):
    line.vars['t_turn_s'] = i*1e-6
    tw = line.twiss(co_guess={'x':-0.07})
    x = (tw.s+157.08/2)%157.08
    y = tw.x*1e3
    mask = np.argsort(x)
    x = x[mask]
    y = y[mask]
    ax.plot(x, y, c=c)
ax.set_xlim(62,98)

norm = Normalize(vmin=min(toplot), vmax=max(toplot))
sm = ScalarMappable(norm = norm, cmap=cm.rainbow)
sm.set_array([])
cbar = ColorbarBase(ax.figure.add_axes([0.99,0.15,0.03,0.85]), cmap=cm.rainbow, norm=norm, orientation='vertical')
cbar.set_label('Turn', fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

f.tight_layout()