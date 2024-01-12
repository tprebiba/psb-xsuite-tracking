#%%
import xtrack as xt
import xobjects as xo
import xfields as xf
import xpart as xp
import numpy as np
import json
import time
import os
from lib.statisticalEmittance import StatisticalEmittance as stE


#%%
#########################################
# Load parameters
#########################################
from simulation_parameters import parameters as p
source_dir = os.getcwd() + '/'
#source_dir = '/afs/cern.ch/user/t/tprebiba/workspace/000_PSB_half-integer_dynamic_crossing_PIC_reference/'
if p['prepare_tune_ramp']:
    with open(source_dir+'time_tables/tunes.json', 'r') as fid:
        d = json.load(fid)
    for key in d:
        p[key] = d[key]


#%%
#########################################
# Load PSB line
#########################################
context = p['context']
line = xt.Line.from_json(source_dir+'psb/psb_line_thin.json')
Cpsb = line.get_length() # 157.08 m
print('Loaded PSB line from psb/psb_line_thin.json.')


#%%
#########################################
# Install space charge nodes
#########################################
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
                    #sigma_z=p['sigma_z'],
                    delta_rms=1e-3,
                    num_spacecharge_interactions=p['num_spacecharge_interactions'])
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
            solver=p['pic_solver'],
            )
    else:
        raise ValueError(f'Invalid mode: {mode}')
    print('Space charge installed')
else:
     print('Skipping space charge...')


#%%
#########################################
# Build tracker
#########################################
line.build_tracker(_context=context)
print('Tracker built')
#line_sc_off = line.filter_elements(exclude_types_starting_with='SpaceCh') # to remove space charge
#print('Keeping line_sc_off: line without space charge knobs.')


#%%
#########################################
# Setup particles for injection
#########################################
print('%s particle distribution.'%p['particle_distribution'])
with open(source_dir+'input/particles_initial.json', 'r') as fid:
    particles_for_injection = xp.Particles.from_dict(json.load(fid), _context=context)
print('Loaded particles from input/particles_initial.json.')
if p['num_injections']==1:
    print('Number of injections = 1.')
    if p['particle_distribution'] == 'simulated':
        particles = particles_for_injection
    elif p['particle_distribution'] == 'real':
        # to be reviewed; current implementation is not correct
        # 'real' particles_for_injection is of length 60000, different from p['n_part']
        particles = particles_for_injection
elif p['num_injections']>1:
    print('Number of injections = %i.'%p['num_injections'])
    
    # Build and insert multi-turn injection element 
    print('Building and inserting multi-turn injection element to PSB lattice.')
    print('Number of injections: ', p['num_injections'])
    print('Number of macroparticles per injection: ', int(p['n_part']/p['num_injections']))
    p_injection = xt.ParticlesInjectionSample(particles_to_inject=particles_for_injection,
                                              line=line,
                                              element_name='injection',
                                              num_particles_to_inject=int(p['n_part']/p['num_injections']))
    line.discard_tracker()
    line.insert_element(index='bi1.tstr1l1', element=p_injection, name='injection')
    line.build_tracker()

    # Generate particle object with unallocated space
    print('Generating particle object with unallocated space.')
    particles = line.build_particles(_capacity=p['n_part']+1, x=0)
    particles.state[0] = -500 # kill the particle added by default
    

#%%
#########################################
# Include injection foil
# I need a to_dict() method for Foil in 
# order to save the line to json
# to be implemented...
#########################################
if p['install_injection_foil']==True:    
    from lib.foil import Foil

    print('Creating PSB foil...')
    thickness = 200 #ug/cm^2
    xmin = -0.099
    #xmax = -0.067 # actual location of the foil
    xmax = +0.099 # for testing
    ymin = -0.029
    ymax = 0.029
    psbfoil = Foil(xmin, xmax, ymin, ymax, thickness)
    print('PSB foil created.')
    print('Setting scatter choice to %s (1: simple (no losses) 0: full (with losses))'%(p['scatterchoice']))
    psbfoil.setScatterChoice(p['scatterchoice'])
    psbfoil.setActivateFoil(1) # activates foil

    print('Inserting foil element to line.')
    line.discard_tracker()
    line.insert_element(index='bi1.tstr1l1', element=psbfoil, name='psbfoil')
    line.build_tracker()


#########################################
# Start lattice at desired location
#########################################
element_to_cycle = p['element_to_cycle']
if element_to_cycle != None:
    line.discard_tracker() # We need to discard the tracker to edit the line
    line.cycle(name_first_element = element_to_cycle, inplace=True)
    print('Changed line starting point to %s.'%(element_to_cycle))
    line.build_tracker()
else:
    print('Line starting point not changed.')

#%%
#########################################
# Last configs
#########################################
if p['prepare_acceleration'] == 0:
    # Switching all RF cavities OFF
    print('Switching off all cavities.')
    line_table = line.get_table()
    element_mask = np.where(line_table['element_type']=='Cavity')[0]
    for i in element_mask:
        cav_name = line_table['name'][i]
        line.element_refs[cav_name].voltage = 0
        line.element_refs[cav_name].frequency = 0
    print('All cavities switched off.')
    line.twiss_default['method'] = '4d'
    print('Twiss method set to 4d.')
line.enable_time_dependent_vars = True
#line.dt_update_time_dependent_vars = 3e-6 # approximately every 3 turns
line.vars.cache_active = False
line.vars['t_turn_s'] = 0.0
output = []
if p['GPU_FLAG']:
    r = stE(context='GPU')
else:
    r = stE(context='CPU')
output=[]
intensity = []


#%%
#########################################
# Start tracking
#########################################
num_turns = p['num_turns']
print('Now start tracking...')
start = time.time()
for ii in range(num_turns):
    print(f'Turn {ii} out of {num_turns}')

    # multi-turn injection + foil
    if p['num_injections']>1:
        if ii == p['num_injections']:
            p_injection.num_particles_to_inject = 0
            print('Injection finished.')
            if p['install_injection_foil']==True:
                #line.psbfoil.setActivateFoil(0) # deactivates foil
                psbfoil.setActivateFoil(0)
                print('Foil deactivated.')
        elif ii<p['num_injections']:
            print('Injecting %i macroparticles.'%(int(p['n_part']/p['num_injections'])))
        intensity.append(particles.weight[particles.state>0].sum())

    # keep particles within the ring circumference
    particles.zeta = (particles.zeta+Cpsb/2)%Cpsb-Cpsb/2

    # track one turn
    #line.track(particles, turn_by_turn_monitor=True)
    line.track(particles, num_turns=1)

    # update output
    bunch_moments=r.measure_bunch_moments(particles)
    output.append([len(r.coordinate_matrix[0]),bunch_moments['nemitt_x'].tolist(),bunch_moments['nemitt_y'].tolist(),bunch_moments['emitt_z'].tolist()])

    # save every some turns
    if ii in p['turns2saveparticles']:
        print(f'Saving turn {ii}')
        with open(source_dir+f'output/particles_turn_{ii:05d}.json', 'w') as fid:
            json.dump(particles.to_dict(), fid, cls=xo.JEncoder)
            print(f'Particles saved to output/particles_turn_{ii:05d}.json.')
        np.save(source_dir+'output/distribution_'+str(int(ii)), r.coordinate_matrix)
        print(f'Distribution saved to output/distribution_{ii}.npy.')
        ouput=np.array(output)
        np.save(source_dir+'output/emittances', output)
        print(f'Emittances saved to output/emittances.npy.')

end = time.time()
print('Tracking finished.')
print('Total seconds = ', end - start)
np.save(source_dir+'output/emittances', output)
print(f'Emittances saved to output/emittances.npy.')
#print('Final x and y emittances: ', bunch_moments['nemitt_x'].tolist()[-1], bunch_moments['nemitt_y'].tolist()[-1])
# %%