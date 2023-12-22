#%%
import xtrack as xt
import xpart as xp
import xobjects as xo
import xfields as xf
import numpy as np
import json
import time
import os

from lib.parabolic_longitudinal_distribution import parabolic_longitudinal_distribution
#from statisticalEmittance.statisticalEmittance import statisticalEmittance as stE
from lib.statisticalEmittance import StatisticalEmittance as stE

GPU_FLAG = False
if GPU_FLAG:
    import cupy as cp


#%%
#############################
# Load parameters
#############################
source_dir = os.getcwd() + '/'
#source_dir = '/afs/cern.ch/user/t/tprebiba/workspace/000_PSB_half-integer_dynamic_crossing_PIC_reference/'
from simulation_parameters import parameters as p
with open(source_dir+'tables/tunes.json', 'r') as fid:
     d = json.load(fid)
for key in d:
     p[key] = d[key]


#%%
#############################
# Load PSB line
#############################
if GPU_FLAG:
    context = xo.ContextCupy()
else:
    context = xo.ContextCpu()
line = xt.Line.from_json(source_dir+'psb/psb_line.json')
Cpsb = line.get_length() # 157.08 m
print('Loaded PSB line')


#%%
#############################
# Install space charge nodes
#############################
space_charge = True
if space_charge:
    mode = 'pic' # 'frozen' or 'pic' or 'quasi-frozen'
    print(f'Installing space charge in {mode} mode')
    # install nodes in lattice frozen-like (exact parameters do not matter if pic is used)
    lprofile = xf.LongitudinalProfileQGaussian(number_of_particles=p['bunch_intensity'], 
                                            sigma_z=p['sigma_z'], 
                                            z0=0, q_parameter=1.)
    xf.install_spacecharge_frozen(line=line,
                    particle_ref=line.particle_ref,
                    longitudinal_profile=lprofile,
                    nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'],
                    sigma_z=p['sigma_z'],
                    num_spacecharge_interactions=p['num_spacecharge_interactions'],
                    tol_spacecharge_position=p['tol_spacecharge_position'])
    if mode == 'frozen':
        pass # Already configured in line
    # switch to pic or quasi-frozen
    elif mode == 'quasi-frozen':
        xf.replace_spacecharge_with_quasi_frozen(line,
                                        update_mean_x_on_track=True,
                                        update_mean_y_on_track=True)
    elif mode == 'pic':
        pic_collection, all_pics = xf.replace_spacecharge_with_PIC(
            _context=context, line=line,
            n_sigmas_range_pic_x=8, # to be reviewed
            n_sigmas_range_pic_y=8, # to be reviewed
            nx_grid=128, ny_grid=128, nz_grid=64, # to be reviewed
            n_lims_x=7, n_lims_y=5,#3,
            #z_range=(-3*p['sigma_z'], 3*p['sigma_z']), 
            z_range=(-Cpsb/2, Cpsb/2), 
            _average_transverse_distribution=True)
    else:
        raise ValueError(f'Invalid mode: {mode}')
    print('Space charge installed')
else:
     print('Skipping space charge...')


#%%
#############################
# Build tracker
#############################
line.build_tracker(_context=context)
print('Tracker built')
line_sc_off = line.filter_elements(exclude_types_starting_with='SpaceCh')
tw = line_sc_off.twiss(particle_ref=line.particle_ref,  at_elements=[0])



#%%
#############################
# Correct chroma
#############################
line.match(vary=[xt.Vary('kbr1xnoh0', step=1e-8)],
           targets = [xt.Target('dqy', 0.0, tol=1e-5)])
print('Chroma corrected')



#%%
#############################
# Generate particles 
#############################
n_part = int(5e5) #int(2e6)
particles = parabolic_longitudinal_distribution(_context=context, num_particles=n_part,
                            total_intensity_particles=p['bunch_intensity'],
                            nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'], sigma_z=p['sigma_z'],
                            particle_ref=line.particle_ref,
                            line=line,
                            #line=line_sc_off
                            )

# "Force" coasting beam
if GPU_FLAG:
    particles.zeta = cp.random.uniform(-Cpsb/2, Cpsb/2, n_part)
else:
    particles.zeta = np.random.uniform(-Cpsb/2, Cpsb/2, n_part)
#particles.delta = np.random.uniform(-1.36e-3, 1.36e-3, n_part)
#particles.delta = np.random.uniform(-1e-10, 1e-10, n_part)
particles0 = particles.copy()
with open(source_dir+'input/particles_initial.json', 'w') as fid:
      json.dump(particles.to_dict(), fid, cls=xo.JEncoder)
print('Particles generated and saved to inputs/particles_initial.json')


#%%
#############################
# Switch OFF RFs
#############################
RF = False
if not RF:
    #line_table = line.get_table()
    #element_mask = np.where(line_table['element_type']=='Cavity')[0]
    #for i in element_mask:
    #    cav_name = line_table['name'][i]
    #    line.element_refs[cav_name].voltage = 0
    #    line.element_refs[cav_name].frequency = 0
    cav_name = 'br.c02'
    line.element_refs[cav_name].voltage = 0
    line.element_refs[cav_name].frequency = 0
    print('RFs switched off')
    try:
        line.twiss()
    except:
         print('Twiss failed (as it should).')
else:
    print('RFs kept on.')


#%%
#############################
# Start tracking
#############################
num_turns = 1
start = time.time()
print('Now start tracking for one turn...')
line.track(particles, num_turns=num_turns, turn_by_turn_monitor=False)
print('Tracking done')
end = time.time()
print(f'Tracking took {end-start} s')
particles1 = particles.copy()


#%%
#############################
# Calculate tune
#############################
def get_phase(u,pu,betu,alfu):
    un = u/np.sqrt(betu)
    pun = u*alfu/np.sqrt(betu) + pu*np.sqrt(betu)
    Ju = un**2 + pun**2
    angle = np.angle(un-1j*pun)
    #angle = np.angle(un-1j*pun)%(2*np.pi)
    #angle = np.angle(un-1j*pun)%(np.pi)
    #angle = np.arctan(pun/un)
    return angle


    
if particles1._num_lost_particles == 0: # no particles were lost
    phasex_0 = get_phase(particles0.x, particles0.px, tw['betx'][0], tw['alfx'][0])
    phasex_1 = get_phase(particles1.x, particles1.px, tw['betx'][0], tw['alfx'][0])
    phasey_0 = get_phase(particles0.y, particles0.py, tw['bety'][0], tw['alfy'][0])
    phasey_1 = get_phase(particles1.y, particles1.py, tw['bety'][0], tw['alfy'][0])
    _qx = (phasex_1 - phasex_0)/(2*np.pi)
    _qy = (phasey_1 - phasey_0)/(2*np.pi)
    #print(_qx)
np.save('output/qx_pic',_qx)
np.save('output/qy_pic',_qy)


# %%
# import matplotlib.pyplot as plt
# plt.close()
# plt.cla()
# #plt.plot(phasex_0, phasey_0, '.', ms=0.03)
# #plt.plot(phasex_1, phasey_1, '.', ms=0.03)
# qxnew = np.copy(_qx)
# qynew = np.copy(_qy)
# qxnew[np.where(qxnew<0)] += 1
# qynew[np.where(qynew<0)] += 1
# plt.plot(_qx, _qy, '.', ms=0.005)
# plt.plot(qxnew, qynew, '.', ms=0.005)
# plt.ylim(0.5,0.6)
# plt.xlim(0,0.2)
# plt.hist2d(_qx, _qy,bins=300, cmap = plt.cm.jet, vmin=10)
# %%
