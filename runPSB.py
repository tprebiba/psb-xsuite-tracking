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

GPU_FLAG = True
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
# Prepare tune ramp
#############################
if p['ramp_tunes_flag']:
     #line.vars['kbrqf'] = line.functions.kbrqf_func(line.vars['t_turn_s'])
     #line.vars['kbrqd'] = line.functions.kbrqd_func(line.vars['t_turn_s'])
     #print('Tune will be ramped.')
     pass
     


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
deactivate_spacecharge = False
if deactivate_spacecharge:
      line = line.filter_elements(exclude_types_starting_with='SpaceCh') # to remove space charge
      print('Space charge deactivated (filtered line)')


#%%
#############################
# Correct chroma
#############################
line.match(vary=[xt.Vary('kbr1xnoh0', step=1e-8)],
           targets = [xt.Target('dqy', 0.0, tol=1e-5)])
print('Chroma corrected')


#%%
#############################
# Add half-integer excitation
#############################
line.vars['kbr1qno816l3'] = p['kbr1qno816l3']
line.vars['kbr1qno412l3'] = p['kbr1qno412l3']
print('Half-integer excitation added')
print('kbr1qno816l3=', p['kbr1qno816l3'])
print('kbr1qno412l3=', p['kbr1qno412l3'])


#%%
#############################
# Generate particles 
#############################
#particles = xp.generate_matched_gaussian_bunch(_context=context, num_particles=p['n_part'],
#                            total_intensity_particles=p['bunch_intensity'],
#                            nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'], sigma_z=p['sigma_z'],
#                            particle_ref=line.particle_ref,
#                            tracker=tracker_sc_off)
particles = parabolic_longitudinal_distribution(_context=context, num_particles=p['n_part'],
                            total_intensity_particles=p['bunch_intensity'],
                            nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'], sigma_z=p['sigma_z'],
                            particle_ref=line.particle_ref,
                            line=line,
                            #line=line_sc_off
                            )
# "Force" coasting beam
if GPU_FLAG:
    particles.zeta = cp.random.uniform(-Cpsb/2, Cpsb/2, p['n_part'])
    print('Forcing coasting beam using cupy.')
else:
    particles.zeta = np.random.uniform(-Cpsb/2, Cpsb/2, p['n_part'])
    print('Forcing coasting beam using numpy.')
#particles.delta = np.random.uniform(-1.36e-3, 1.36e-3, n_part) # not parabolic
with open(source_dir+'input/particles_initial.json', 'w') as fid:
      json.dump(particles.to_dict(), fid, cls=xo.JEncoder)
print('Number of macroparticles: ', p['n_part'])
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
# Last configs
#############################
#line.enable_time_dependent_vars = True
#line.dt_update_time_dependent_vars = 3e-6 # approximately every 3 turns
#line.vars.cache_active = False
output = []
if GPU_FLAG:
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
num_turns = p['num_turns']
print('Now start tracking...')
start = time.time()
for ii in range(num_turns):
      print(f'Turn {ii} out of {num_turns}')

      # ramp quadrupoles: not needed when using the kbrqf_func and kbrqd_func
      line.vars['kbrqf'] = p['kf_ini'] + (p['kf_fin']-p['kf_ini'])*ii/num_turns
      line.vars['kbrqd'] = p['kd_ini'] + (p['kd_fin']-p['kd_ini'])*ii/num_turns
      print('Model Qx = ', p['qx_ini'] + (p['qx_fin']-p['qx_ini'])*ii/num_turns)
      print('Model Qy = ', p['qy_ini'] + (p['qy_fin']-p['qy_ini'])*ii/num_turns)

      # keep particles within circumference
      particles.zeta = (particles.zeta+Cpsb/2)%Cpsb-Cpsb/2
      
      # track one turn
      line.track(particles, turn_by_turn_monitor=True)
      #t = line.twiss(method='4d')
      #print(t.qx, t.qy)
      
      # update output
      bunch_moments=r.measure_bunch_moments(particles)
      output.append([len(r.coordinate_matrix[0]),bunch_moments['nemitt_x'].tolist(),bunch_moments['nemitt_y'].tolist(),bunch_moments['emitt_z'].tolist()])
      
      # save every some turns
      ii2save = [100, 200, 1000, 10000, 20000, 30000, 40000, 50000]
      #ii2save = [10000, 14000, 16000, 18000, 20000, 22000, 24000, 26000, 28000, 30000, 32000, 34000]
      #ii2save = [1000, 200, 300, 40, 50, 100, 2000, 3000]
      if (ii in ii2save) or (ii > 50000):
            if ii%1000==0:
            #if True:
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
