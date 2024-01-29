import numpy as np
import matplotlib.pyplot as plt


def plot_phasespace(particles, turn, png_dir, 
                    bins=500, vmin=0.1,
                    xy_extent_mm=22, pxpy_extent_mrad=7, z_extent_m=157.08/2, dE_extent_MeV=2,
                    fontsize=10,
                    GPU_FLAG=False, mask_lost_particles=True):

    x = particles.x
    y = particles.y
    xp = particles.px
    yp = particles.py
    z = particles.zeta
    dE = particles.ptau*particles.p0c/1e6
    state = particles.state
    print('number of particles in turn ', turn, 'is ', len(x))

    if GPU_FLAG:
        x = x.get()
        y = y.get()
        xp = xp.get()
        yp = yp.get()
        z = z.get()
        dE = dE.get()
        state = state.get()
    
    if mask_lost_particles:
        mask = particles.state > 0
        x = x[mask]
        y = y[mask]
        xp = xp[mask]
        yp = yp[mask]
        z = z[mask]
        dE = dE[mask]

    my_cmap = plt.cm.jet
    my_cmap.set_under('w',0.1)

    f, axs = plt.subplots(2,2,figsize=(7,5),facecolor='white')
    f.suptitle('Turn %s'%turn, fontsize=fontsize)
    
    ax = axs[0,0]
    ax.hist2d(x*1000., xp*1000.,bins=bins, cmap=my_cmap, vmin=vmin)
    #ax.plot(x*1000., xp*1000.,ms=0.5,color='blue')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('xp [mrad]')
    ax.set_xlim(-xy_extent_mm,xy_extent_mm)
    ax.set_ylim(-pxpy_extent_mrad,pxpy_extent_mrad)
    
    ax = axs[0,1]
    ax.hist2d(y*1000., yp*1000.,bins=bins, cmap=my_cmap, vmin=vmin)
    #ax.plot(y*1000., yp*1000.,ms=0.5,color='blue')
    ax.set_xlabel('y [mm]')
    ax.set_ylabel('yp [mrad]')
    ax.set_xlim(-xy_extent_mm,xy_extent_mm)
    ax.set_ylim(-pxpy_extent_mrad,pxpy_extent_mrad)
    
    ax = axs[1,0]
    ax.hist2d(x*1000., y*1000.,bins=bins, cmap=my_cmap, vmin=vmin)
    #ax.plot(x*1000., y*1000.,ms=0.5,color='blue')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm] ')
    ax.set_xlim(-xy_extent_mm,xy_extent_mm)
    ax.set_ylim(-xy_extent_mm,xy_extent_mm)
    
    ax = axs[1,1]
    ax.hist2d(z, dE,bins=bins, cmap=my_cmap, vmin=vmin)
    #ax.plot(z, dE, '.', c='blue')
    ax.set_xlabel('z [m]')
    ax.set_ylabel('dE [MeV] ')
    ax.set_xlim(-z_extent_m,z_extent_m)
    ax.set_ylim(-dE_extent_MeV,dE_extent_MeV)
    
    for ax in axs.flatten():
        ax.xaxis.label.set_size(fontsize)
        ax.yaxis.label.set_size(fontsize)
        ax.tick_params(labelsize=fontsize)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #plt.show()
    plt.savefig('%sphasespace_%s.png'%(png_dir, turn), dpi=100)
    plt.close()