import xtrack as xt
import xpart as xp
import xobjects as xo
from simulation_parameters import parameters as p
from lib.parabolic_longitudinal_distribution import parabolic_longitudinal_distribution
import numpy as np
import json
if p['GPU_FLAG']:
    import cupy as cp

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('007_generate_particle_distribution.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')


#########################################
# Load PSB line in xsuite
#########################################
print('Loading PSB line from psb/psb_line_thin.json.')
line = xt.Line.from_json('psb/psb_line_thin.json') # to get the lengths of the injection bumpers
line.build_tracker()
Cpsb = line.get_length() # 157.08 m


#############################
# Generate particles 
#############################
print('Generating particles...')
context = p['context']
if p['longitudinal_shape'] == 'gaussian':
    particles = xp.generate_matched_gaussian_bunch(_context=context, num_particles=p['n_part'],
                            total_intensity_particles=p['bunch_intensity'],
                            nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'], sigma_z=p['sigma_z'],
                            particle_ref=line.particle_ref,
                            line=line
                            )
    print('Gaussian longitudinal distribution.')
elif p['longitudinal_shape'] == 'parabolic':
    particles = parabolic_longitudinal_distribution(_context=context, num_particles=p['n_part'],
                                total_intensity_particles=p['bunch_intensity'],
                                nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'], sigma_z=p['sigma_z'],
                                particle_ref=line.particle_ref,
                                line=line,
                                #line=line_sc_off
                                )
    print('Parabolic longitudinal distribution.')
elif p['longitudinal_shape'] == 'coasting':
    # to be improved...
    particles = parabolic_longitudinal_distribution(_context=context, num_particles=p['n_part'],
                                total_intensity_particles=p['bunch_intensity'],
                                nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'], sigma_z=p['sigma_z'],
                                particle_ref=line.particle_ref,
                                line=line,
                                #line=line_sc_off
                                )
    # "Force" coasting beam
    if p['GPU_FLAG']:
        particles.zeta = cp.random.uniform(-Cpsb/2, Cpsb/2, p['n_part'])
        print("'Forcing' coasting beam using cupy.")
    else:
        particles.zeta = np.random.uniform(-Cpsb/2, Cpsb/2, p['n_part'])
        print("'Forcing' coasting beam using numpy.")
    print('Coasting beam.')
#particles.delta = np.random.uniform(-1.36e-3, 1.36e-3, n_part) # not parabolic
with open('input/particles_initial.json', 'w') as fid:
      json.dump(particles.to_dict(), fid, cls=xo.JEncoder)
print('Number of macroparticles: ', p['n_part'])
print('Particles generated and saved to inputs/particles_initial.json.')