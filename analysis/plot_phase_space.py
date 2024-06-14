#%%
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xplt
import matplotlib.gridspec as gridspec

#%%
png_dir = './'
#source_dir = '../input/'
source_dir = '../output/'
files = glob.glob(source_dir+'*.json')
#files = files[0:15]
#print(files)


#%%
##############################################
# Plot phase space a la pyorbit
##############################################
for file in files:
#for file in files[10:20]:
    print(file)

    df = pd.read_json(file)

    try:
         turn = int(file.split('_')[-1].split('.')[0])
    except:
         turn = '00'
    x = df['x']
    y = df['y']
    xp = df['px']
    yp = df['py']
    z = df['zeta']
    dE = df['ptau']*df['p0c']/1e6

    mask = np.where(df['state']>0)[0]
    #mask = np.where(df['state']<0)[0]
    #mask = np.where(df['state']!=np.nan)[0]
    x = x[mask]
    y = y[mask]
    xp = xp[mask]
    yp = yp[mask]
    z = z[mask]
    dE = dE[mask]
    print('Number of particles: %s'%len(x))
    print('Number of lost particles: %s'%(len(df['x'])-len(x)))

    bins=200 #500
    my_cmap = plt.cm.jet
    my_cmap.set_under('w',0.1)
    #f, axs = plt.subplots(2,2,figsize=(10,8),facecolor='white')
    #fontsize=15
    f, axs = plt.subplots(2,2,figsize=(7,5),facecolor='white')
    fontsize=10
    f.suptitle('Turn %s'%turn, fontsize=fontsize)
    vmin = 0.1 # the smaller it is, the less sensitive is the colormap
    ax = axs[0,0]
    ax.hist2d(x*1000., xp*1000.,bins=bins, cmap=my_cmap, vmin=vmin)
    #ax.plot(x*1000., xp*1000.,ms=0.5,color='blue')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('xp [mrad]')
    #ax.set_xlim(-22,22)
    ax.set_ylim(-7,7)
    ax.set_xlim(-90,22)
    ax = axs[0,1]
    ax.hist2d(y*1000., yp*1000.,bins=bins, cmap=my_cmap, vmin=vmin)
    #ax.plot(y*1000., yp*1000.,ms=0.5,color='blue')
    ax.set_xlabel('y [mm]')
    ax.set_ylabel('yp [mrad]')
    ax.set_xlim(-22,22)
    ax.set_ylim(-7,7)
    ax = axs[1,0]
    ax.hist2d(x*1000., y*1000.,bins=bins, cmap=my_cmap, vmin=vmin)
    #ax.plot(x*1000., y*1000.,ms=0.5,color='blue')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm] ')
    #ax.set_xlim(-22,22)
    #ax.set_ylim(-22,22)
    ax.set_xlim(-90,22)
    ax.set_ylim(-22,22)
    ax = axs[1,1]
    #ax.hist2d(z, dE,bins=bins, cmap=my_cmap, vmin=vmin)
    ax.plot(z, dE,'.', c='blue', ms=1)
    ax.set_xlabel('z [m]')
    ax.set_ylabel('dE [MeV] ')
    #ax.set_ylim(-2,2)
    #ax.set_xlim(-157.2/2, 157.2/2)
    for ax in axs.flatten():
            ax.xaxis.label.set_size(fontsize)
            ax.yaxis.label.set_size(fontsize)
            ax.tick_params(labelsize=fontsize)
    #plt.suptitle('turn %s'%turn, fontsize=fontsize)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    #plt.savefig('%s/phasespace_%s.png'%(png_dir, turn), dpi=100)
    #plt.close(f)


# %%
##############################################
# Plot phase space a la xplt
##############################################    
for file in files[-3:-1]:
     df = pd.read_json(file)
     plot = xplt.PhaseSpacePlot(df, kind='x')
     plot.fig.tight_layout()
     plot.fig.show()


#%%
##############################################
# Plot phase space a la tirsi
##############################################
for file in files[0:30]:
     df = pd.read_json(file)
     turn = int(file.split('_')[-1].split('.')[0])

     my_cmap = plt.cm.jet
     my_cmap.set_under('w',1)
     fontsize=20
     bins = 200
     vmin = 0.1
     #vmin = 1

     fig = plt.figure(figsize=(8,6), facecolor='white')
     gs = gridspec.GridSpec(6, 6)
     ax = plt.subplot(gs[1:6, 0:4])
     ax_xDist = plt.subplot(gs[0, 0:4],sharex=ax)
     ax_yDist = plt.subplot(gs[1:6, 4],sharey=ax)

     ax.hist2d(df.x*1e3, df.px*1e3,
               bins=500, cmap = my_cmap, vmin=vmin, 
               #range=[[min(df.x)*1e3, max(df.x)*1e3], [min(df.px)*1e3, max(df.px)*1e3]]#, norm=mcolors.PowerNorm(gamma)
               range = [[-90, 22], [-7, 7]]
               )

     #r.plot_resonance(ax)
     ax.set_xlabel('$x$ (mm)', fontsize=fontsize)
     ax.set_ylabel('$p_x$ (mrad)', fontsize=fontsize)
     ax.tick_params(axis='both', labelsize=fontsize)
     ax.set_xlim(-90, 22)
     ax.set_ylim(-7, 7)

     ax_xDist.hist(df.x*1e3,bins=bins,align='mid', color='black', range=(-90, 22))
     ax_xDist.tick_params(axis='x', labelleft=False, labelright=False, labeltop=False, labelbottom=False)
     ax_xDist.tick_params(axis='y', labelleft=False, labelright=False, labeltop=False, labelbottom=False)
     nsigmas = 2.5
     #ax_xDist.axvline(np.mean(qxnew)-nsigmas*np.std(qxnew), color='red', ls='-')
     #ax_xDist.axvline(np.mean(qxnew)+nsigmas*np.std(qxnew), color='red', ls='-')
     #print('Total tune spread in x: ', nsigmas*np.std(qxnew)*2)
     #ax_xDist.set(ylabel='count')
     #ax_xCumDist = ax_xDist.twinx()
     #ax_xCumDist.hist(qxnew,bins=bins,cumulative=True,histtype='step',density=True,color='r',align='mid')
     #ax_xCumDist.tick_params('y', colors='r')
     #ax_xCumDist.set_ylabel('cumulative',color='r')

     ax_yDist.hist(df.px*1e3,bins=bins,orientation='horizontal',align='mid', color='black', range=(-7, 7))
     ax_yDist.tick_params(axis='x', labelleft=False, labelright=False, labeltop=False, labelbottom=False)
     ax_yDist.tick_params(axis='y', labelleft=False, labelright=False, labeltop=False, labelbottom=False)
     nsigmas = 2.5
     #ax_yDist.axhline(np.mean(qynew)-nsigmas*np.std(qynew), color='red', ls='-')
     #ax_yDist.axhline(np.mean(qynew)+nsigmas*np.std(qynew), color='red', ls='-')
     #print('Total tune spread in y: ', nsigmas*np.std(qynew)*2)
     #ax_yDist.set(xlabel='count')
     #ax_yCumDist = ax_yDist.twiny()
     #ax_yCumDist.hist(qynew,bins=bins,cumulative=True,histtype='step',density=True,color='r',align='mid',orientation='horizontal')
     #ax_yCumDist.tick_params('x', colors='r')
     #ax_yCumDist.set_xlabel('cumulative',color='r')

     plt.suptitle('Turn %s'%turn, fontsize=fontsize)
     fig.tight_layout()
     plt.show()
     fig.savefig('%s/phasespace_%s.png'%(png_dir, turn), dpi=100)
# %%
