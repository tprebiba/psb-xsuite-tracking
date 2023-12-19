#%%
import xtrack as xt
import xobjects as xo
import xfields as xf
import xpart as xp
import numpy as np
import json
import time
import os

#from statisticalEmittance.statisticalEmittance import statisticalEmittance as stE
from lib.statisticalEmittance import StatisticalEmittance as stE


#%%
#############################
# Load parameters
#############################
from simulation_parameters import parameters as p
source_dir = os.getcwd() + '/'
#source_dir = '/afs/cern.ch/user/t/tprebiba/workspace/000_PSB_half-integer_dynamic_crossing_PIC_reference/'
if p['prepare_tune_ramp']:
    with open(source_dir+'time_tables/tunes.json', 'r') as fid:
        d = json.load(fid)
    for key in d:
        p[key] = d[key]


#%%
#############################
# Load PSB line
#############################
context = p['context']
line = xt.Line.from_json(source_dir+'psb/psb_line_thin.json')
Cpsb = line.get_length() # 157.08 m
print('Loaded PSB line from psb/psb_line_thin.json.')


#%%
#############################
# Install space charge nodes
#############################
if p['install_space_charge']:
    mode = p['space_charge_mode']
    print(f'Installing space charge in {mode} mode')
    # install nodes in lattice frozen-like (exact parameters do not matter if pic is used)
    lprofile = xf.LongitudinalProfileQGaussian(number_of_particles=p['bunch_intensity'], 
                                            sigma_z=p['sigma_z']*10/2, 
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
            _average_transverse_distribution=False)
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
line_sc_off = line.filter_elements(exclude_types_starting_with='SpaceCh') # to remove space charge
print('Keeping line_sc_off: line without space charge knobs.')


#%%
#############################
# Get particles from input
#############################
if p['particle_distribution'] == 'simulated':
    print('Simulated particle distribution.')
    with open(source_dir+'input/particles_initial.json', 'r') as fid:
        particles = xp.Particles.from_dict(json.load(fid), _context=context)
    print('Loaded particles from input/particles_initial.json.')
elif p['particle_distribution'] == 'real':
    print('Real particle distribution.')
    with open('input/part_for_injection.json', 'r') as fid:
        part_for_injection = xp.Particles.from_dict(json.load(fid), _context=context)
    print('Loaded particles from input/part_for_injection.json.')


#%%
#############################
# Multi-turn injection
#############################
if p['particle_distribution'] == 'real':
    ###############################
    # Build and insert multi-turn
    # injection element 
    ###############################
    print('Building and inserting multi-turn injection element to PSB lattice.')
    print('Number of injections: ', p['num_injections'])
    print('Number of particles per injection: ', p['n_part']/p['num_injections'])
    p_injection = xt.ParticlesInjectionSample(particles_to_inject=part_for_injection,
                                              line=line,
                                              element_name='injection',
                                              num_particles_to_inject=p['n_part']/p['num_injections'])
    line.discard_tracker()
    line.insert_element(index='bi1.tstr1l1', element=p_injection, name='injection')
    line.build_tracker()

    ###############################
    # Generate particle object with 
    # unallocated space
    ###############################
    particles = line.build_particles(_capacity=p['n_part'], x=0)
    particles.state[0] = -500 # kill the particle added by default


#%%
#############################
# Last configs
#############################
line.enable_time_dependent_vars = True
#line.dt_update_time_dependent_vars = 3e-6 # approximately every 3 turns
line.vars.cache_active = False
line.vars['t_turn_s'] = 0.0
output = []
if p['GPU_FLAG']:
    r = stE(context='GPU')
else:
    r = stE(context='CPU')
bunch_moments=r.measure_bunch_moments(particles)
print(bunch_moments['nemitt_x'])
print(bunch_moments['nemitt_y'])
output=[]
#np.save(source_dir+'output/distribution_thin_seq', r.coordinate_matrix)


#%%
#############################
# Start tracking
#############################
intensity = []
num_turns = p['num_turns']
print('Now start tracking...')
start = time.time()
for ii in range(num_turns):
      print(f'Turn {ii} out of {num_turns}')

      # multi-turn injection
      if p['particle_distribution'] == 'real':
        if ii == p['num_injections']:
            p_injection.num_particles_to_inject = 0
        print('Stopping...')
        intensity.append(particles.weight[particles.state>0].sum())

      # keep particles within circumference
      #particles.zeta = (particles.zeta+Cpsb/2)%Cpsb-Cpsb/2
      
      # track one turn
      line.track(particles, turn_by_turn_monitor=True)
      #t = line.twiss(method='4d')
      #print(t.qx, t.qy)
      
      # update output
      bunch_moments=r.measure_bunch_moments(particles)
      output.append([len(r.coordinate_matrix[0]),bunch_moments['nemitt_x'].tolist(),bunch_moments['nemitt_y'].tolist(),bunch_moments['emitt_z'].tolist()])
      
      # save every some turns
      if ii in p['turns2saveparticles']:
        print(f'save turn {ii}')
        with open(source_dir+f'output/particles_turn_{ii:05d}.json', 'w') as fid:
            json.dump(particles.to_dict(), fid, cls=xo.JEncoder)
        np.save(source_dir+'output/distribution_'+str(int(ii)), r.coordinate_matrix)
        ouput=np.array(output)
        np.save(source_dir+'output/emittances', output)

end = time.time()
bunch_moments=r.measure_bunch_moments(particles)
print('epsn_x = ',bunch_moments['nemitt_x'])
print('epsn_y = ',bunch_moments['nemitt_y'])
print('eps_z = ',bunch_moments['emitt_z'])
print('time = ', end - start)
ouput=np.array(output)
np.save(source_dir+'output/emittances', output)

# %%
