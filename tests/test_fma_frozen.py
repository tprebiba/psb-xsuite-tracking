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
print('Loaded PSB line')   


#%%
#############################
# Install space charge nodes
#############################
space_charge = True
if space_charge:
    mode = 'frozen' # 'frozen' or 'pic' or 'quasi-frozen'
    print(f'Installing space charge in {mode} mode')
    # install nodes in lattice frozen-like (exact parameters do not matter if pic is used)
    lprofile = xf.LongitudinalProfileQGaussian(number_of_particles=p['bunch_intensity'], 
                                           sigma_z=p['sigma_z']*10/3, 
                                           z0=0, q_parameter=1.)
    #lprofile = xf.LongitudinalProfileCoasting(context=context,beam_line_density=0)
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
            n_lims_x=7, n_lims_y=3,
            z_range=(-3*p['sigma_z'], 3*p['sigma_z']), 
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
deactivate_spacecharge = False
if deactivate_spacecharge:
      line = line.filter_elements(exclude_types_starting_with='SpaceCh') # to remove space charge
      print('Space charge deactivated (filtered line)')


#%%
#############################
# Correct chroma
#############################
#line.match(vary=[xt.Vary('kbr1xnoh0', step=1e-8)],
#           targets = [xt.Target('dqy', 0.0, tol=1e-5)])
#print('Chroma corrected')



#%%
#############################
# Generate particles 
#############################
nrthetas = 20#-10
nrrhos = 30#-10
x_norm, y_norm, _, _ = xp.generate_2D_polar_grid(
    theta_range=(0.01, np.pi/2-0.01),
    ntheta = nrthetas,
    r_range = (0.1, 7),
    nr = nrrhos)

particles = xp.build_particles(line=line, particle_ref=line.particle_ref,
                               x_norm=x_norm, y_norm=y_norm, delta=0,
                               nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'])
npart = particles._capacity

# "Force" coasting beam
Cpsb = line.get_length() # 157.08 m
#if GPU_FLAG:
#    particles.zeta = cp.random.uniform(-Cpsb/2, Cpsb/2, p['n_part'])
#else:
    #particles.zeta = np.random.uniform(-Cpsb/2, Cpsb/2, p['n_part'])
#    particles.zeta = np.random.uniform(-Cpsb/2, Cpsb/2, 20*30)
#particles.delta = np.random.uniform(-1.36e-3, 1.36e-3, 20*30)
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
num_turns = 1200 #p['num_turns']
print('Now start tracking...')
#lasttrack_x = np.zeros((npart, num_turns))
#lasttrack_y = np.zeros((npart, num_turns))
#for i in range(num_turns):
#    print(i)
#    line.track(particles, num_turns=1, turn_by_turn_monitor=False)
#    lasttrack_x[:,i] = particles.x.copy()
#    lasttrack_y[:,i] = particles.y.copy()
line.track(particles, num_turns=num_turns, turn_by_turn_monitor=True)
np.save('output/x',line.record_last_track.x)
#np.save('output/x',lasttrack_x)
np.save('output/y',line.record_last_track.y)
#np.save('output/y',lasttrack_y)
#np.save('output/px',line.record_last_track.px)
#np.save('output/py',line.record_last_track.py)
np.save('output/z',line.record_last_track.zeta)
np.save('output/dp',line.record_last_track.delta)
print("Done!")
# %%
