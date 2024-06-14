import xtrack as xt
import xpart as xp
import xobjects as xo
from simulation_parameters import parameters as p
from lib.parabolic_longitudinal_distribution import parabolic_longitudinal_distribution
import numpy as np
import json
import pandas as pd

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('007_generate_particle_distribution.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

#########################################
# Load PSB line in xsuite and cycle at foil
#########################################
print('Loading PSB line from psb/psb_line_thin.json.')
line = xt.Line.from_json('psb/psb_line_thin.json')
line.build_tracker()
Cpsb = line.get_length() # 157.08 m

#########################################
# Generate particle distribution
#########################################
if p['particle_distribution'] == 'simulated':
    print('Simulated particle distribution.')

    element_to_cycle = p['element_to_cycle']
    if element_to_cycle != None:
        print('Element to cycle: %s.'%element_to_cycle)
        print('Changing starting point of line to generate (simulated) matched distribution with correct optics.')
        line.discard_tracker() # We need to discard the tracker to edit the line
        line.cycle(name_first_element = p['element_to_cycle'], inplace=True)
        # If injection perturbations are included, and line is cycled to the foil (bi1.tstr1l1), line.twiss() needs
        # to be called either with co_search_at some other location (for example 'psb1$start') or to give a good 
        # initial guess of the closed orbit at the foil. This is because of the large bump at the foil (~81 mm).
        # Thus, the particle generators below will fail...
        line.build_tracker()
        line.to_json('psb/psb_line_thin.json')
        print('Line saved to psb/psb_line_thin.json')

    print('Generating particles...')
    context = xo.ContextCpu()
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
                                    line=line
                                    )
        print('Parabolic longitudinal distribution.')
    elif p['longitudinal_shape'] == 'coasting':
        # To be improved...
        # Generate particles
        particles = parabolic_longitudinal_distribution(_context=context, num_particles=p['n_part'],
                                    total_intensity_particles=p['bunch_intensity'],
                                    nemitt_x=p['nemitt_x'], nemitt_y=p['nemitt_y'], sigma_z=p['sigma_z'],
                                    particle_ref=line.particle_ref,
                                    line=line,
                                    #line=line_sc_off
                                    )
        
        # "Force" coasting beam
        particles.zeta = np.random.uniform(-Cpsb/2, Cpsb/2, p['n_part'])
        #particles.delta = np.random.uniform(-1.36e-3, 1.36e-3, n_part) # not parabolic
        print("'Forcing' coasting beam using numpy (uniform zeta).")

    # Add injection missteering
    particles.x += p['injection_missteering_x']
    particles.y += p['injection_missteering_y']
    
    with open('input/particles_initial.json', 'w') as fid:
        json.dump(particles.to_dict(), fid, cls=xo.JEncoder)
    print('Number of macroparticles: ', p['n_part'])
    print('Particles generated and saved to inputs/particles_initial.json.')

elif p['particle_distribution'] == 'real':
    print('Real particle distribution.')
    
    print('Twissing.')
    tw = line.twiss()
    #tw = line.twiss(co_search_at='psb1$start') # closed orbit at PSB start, otherwise need to give co_guess (xtrack 0.52.0)
    co_x_at_foil = tw['x', 'bi1.tstr1l1_entry']
    co_y_at_foil = tw['y', 'bi1.tstr1l1_entry']
    print('Closed orbit at foil: x = %s m, y = %s m.'%(co_x_at_foil, co_y_at_foil))

    # Multi-turn injection is performed at the foil; mismatched otherwise
    line.discard_tracker() # We need to discard the tracker to edit the line
    line.cycle(name_first_element = 'bi1.tstr1l1', inplace=True)
    line.build_tracker()
    print('Changed line starting point to bi1.tstr1l1 (foil).')
    line.to_json('psb/psb_line_thin.json')
    print('Line saved to psb/psb_line_thin.json')
    
    fname = 'L4_particle_distribution/atPSBfoil-450keV.txt'
    #fname = 'L4_particle_distribution/250keV.txt'
    print('Reading %s.'%fname)
    df = pd.read_table(fname, skiprows=3,
        names="x x' y y' z z' Phase Time Energy Loss".split())
    
    # read first and second row of file
    with open(fname) as fid:
        line1 = fid.readline()
        line2 = fid.readline()
    keys = line1.split()
    values = line2.split()
    linac4parameters = dict(zip(keys, values))
    
    kin_energy_ev = df.Energy.values * 1e6
    tot_energy_ev = kin_energy_ev + xt.PROTON_MASS_EV
    p0c = line.particle_ref.p0c[0]
    #tot_energy0_ev = line.particle_ref.energy0[0] # taking energy from reference particle: gives slight energy mismatch
    tot_energy0_ev = float(linac4parameters['Beam_energy(MeV)'])*1e6+line.particle_ref.mass0 # taking energy from Linac4 file
    ptau = (tot_energy_ev - tot_energy0_ev) / p0c
    print('Building particles object from %s.'%fname)
    particles = xt.Particles(q0=1, mass0=xt.PROTON_MASS_EV, p0c=line.particle_ref.p0c[0],ptau=ptau)
    particles.x = df.x.values * 1e-3
    particles.y = df.y.values * 1e-3
    
    # injection at closed orbit
    particles.x += co_x_at_foil
    particles.y += co_y_at_foil

    # add injection missteering
    particles.x += p['injection_missteering_x']
    particles.y += p['injection_missteering_y']

    # Following analysis from PyOrbit
    # TO BE REVIEWED
    z_lim = p['choppingFactor']*2*np.pi*25
    z_offset = p['z_offset']
    z_min, z_max,  = -z_lim/2.+z_offset, z_lim/2.+z_offset 
    speed_of_light = 299792458 # m/s
    Linac4_bucket_length = line.particle_ref._beta0[0]*speed_of_light/(1e6*float(linac4parameters['Beam_Frequency(MHz)']))
    bucket_min = np.round((z_min+Linac4_bucket_length/2.)/Linac4_bucket_length)
    bucket_max = np.round((z_max-Linac4_bucket_length/2.)/Linac4_bucket_length)
    particles.zeta = df.z.values * 1e-3 + np.int_(np.round(0.5+np.random.uniform(bucket_min-1,bucket_max,len(df.z.values))))*Linac4_bucket_length
    
    particles.px = df["x'"].values * 1e-3 * (1 + particles.delta)
    particles.py = df["y'"].values * 1e-3 * (1 + particles.delta)
    particles.weight = p['macrosize']
    print('Macrosize = ', p['macrosize'])
    with open('input/particles_initial.json', 'w') as fid:
        json.dump(particles.to_dict(), fid, cls=xo.JEncoder)
    print('Saved %s to particles_initial.json.'%fname)