import xobjects as xo
import numpy as np

#######################################
# 1) Tracking parameters
#######################################
#num_turns = 100000 # number of turns to track
num_turns = 500
turns2saveparticles = [1,2,5,10,20,50,100,200, 300, 500, 1000, 2000, num_turns-1]

#######################################
# 2) Beam parameters
#######################################
#n_part = int(5e5) #int(5e5) # number of macroparticles: 2e5 were used in pyorbit
n_part = int(5e3)
bunch_intensity = 40.0e10 # number of particles per bunch
macrosize = bunch_intensity/n_part # number of particles (charges) per macroparticle

particle_distribution = 'real' # 'simulated' or 'real
if particle_distribution == 'simulated':
    nemitt_x = 1.0e-6 # normalized emittance in x
    nemitt_y = 1.0e-6 # normalized emittance in y
    sigma_z = (400/4)*0.525*0.3 # bunch length in m
    longitudinal_shape = 'parabolic' # 'parabolic' or 'coasting' or 'gaussian'
elif particle_distribution == 'real':
    nemitt_x = 0.4e-6 # guess, to be revised
    nemitt_y = 0.3e-6 # guess, to be revised
    sigma_z = 10 # guess, to be revised
    longitudinal_shape = np.nan # not used

num_injections = 1 # if > 1: a multi-turn injection is setup

qx_ini = 4.40
qy_ini = 4.45

#######################################
# 3) Lattice imperfections
#######################################
correct_chroma = False # if True, correct chromaticity
chroma_plane = 'y' # to correct horizontal ('x') or vertical ('y') chromaticity

# Half-integer excitation (deltaI_816 = -2A)
include_field_errors = False
field_errors = {
    'kbr1qno816l3': -6.15363e-4,#*10 # half-integer excitation
    'kbr1qno412l3': 0,
}

#######################################
# 4) Flags to control simulation flow
#######################################
include_injection_chicane = 0 # if 1, 002A_include_injection_chicane.py is executed
on_chicane_k0 = 1 # if 1, activates edge focusing of injection chicane
on_chicane_k2 = 1 # if 1, activates eddy currents of injection chicane

include_injection_chicane_correction = 0 # if 1, 002B_include_injection_chicane_correction.py is executed
on_chicane_tune_corr = 1 # if 1, activates tune correction of injection chicane
on_chicane_beta_corr = 1 # if 1, activates beta correction of injection chicane

prepare_acceleration = 2 # 0: all RF OFF, 1: nominal PSB acceleration (double RF), 2: flat bottom (single RF at 8kV)
twiss_mode = '4d' # '4d' or '6d
zeta0 = 17.5 # if double RF, 6d-twiss at zeta0=0 fails because is unstable fixed point; this is a guess of the stable fixed point

slices = 3 # number of slices in thin lattice
# to have the starting point of the lattice at a different location, None otherwise
element_to_cycle = None # line starts at the start of the 1st sector (NOT at the stripping foil)
#element_to_cycle = 'bi1.tstr1l1' # stripping foil
#element_to_cycle = 'br1.bwsv11l1' # vertical LIU wire scanner

prepare_tune_ramp = 0 # if 1, 004_prepare_tune_ramp.py is executed
on_tune_ramp = 1 # if 1, activates tune ramp
qx_fin = 4.17 # final horizontal tune
qy_fin = 4.23 # final vertical tune

install_space_charge = True # if True, space charge is installed
space_charge_mode = 'frozen' # 'frozen' or 'pic' or 'quasi-frozen'
num_spacecharge_interactions = 160 # space charge interactions per turn
pic_solver = 'FFTSolver2p5DAveraged' # `FFTSolver2p5DAveraged` or `FFTSolver2p5D`

install_injection_foil = False # if True, injection foil is installed
scatterchoice = 1 # 1: simple (no losses) 0: full (with losses) not working (for now)!

GPU_FLAG = False # if True, GPU is used
if GPU_FLAG:
    context = xo.ContextCupy()
else:
    context = xo.ContextCpu()

#######################################
# Store parameters in dictionary
#######################################
parameters = {
    # Tracking parameters
    'num_turns': num_turns,
    'turns2saveparticles': turns2saveparticles,

    # Beam parameters
    'n_part': n_part,
    'bunch_intensity': bunch_intensity,
    'macrosize': macrosize,
    'particle_distribution': particle_distribution,
    'nemitt_x': nemitt_x,
    'nemitt_y': nemitt_y,
    'sigma_z': sigma_z,
    'longitudinal_shape': longitudinal_shape,
    'num_injections': num_injections,
    'qx_ini': qx_ini,
    'qy_ini': qy_ini,

    # Lattice imperfections
    'correct_chroma': correct_chroma,
    'chroma_plane': chroma_plane,
    'include_field_errors': include_field_errors,
    'field_errors': field_errors, # is a dictionary
    
    # Flags to control simulation flow
    'include_injection_chicane': include_injection_chicane,
    'on_chicane_k0': on_chicane_k0,
    'on_chicane_k2': on_chicane_k2,
    'include_injection_chicane_correction': include_injection_chicane_correction,
    'on_chicane_tune_corr': on_chicane_tune_corr,
    'on_chicane_beta_corr': on_chicane_beta_corr,
    'prepare_acceleration': prepare_acceleration,
    'twiss_mode': twiss_mode,
    'zeta0': zeta0,
    'slices': slices,
    'element_to_cycle': element_to_cycle,
    'prepare_tune_ramp': prepare_tune_ramp,
    'on_tune_ramp': on_tune_ramp,
    'qx_fin': qx_fin,
    'qy_fin': qy_fin,
    'install_space_charge': install_space_charge,
    'space_charge_mode': space_charge_mode,
    'num_spacecharge_interactions': num_spacecharge_interactions,
    'pic_solver': pic_solver,
    'install_injection_foil': install_injection_foil,
    'scatterchoice': scatterchoice,
    'GPU_FLAG': GPU_FLAG,
    'context': context,
}
