import xobjects as xo
import numpy as np

p = {} 

# Tracking parameters
p['num_turns'] = 2000 # number of turns to track
p['turns2saveparticles'] = [1,2,3,4,5,6,7,8,9,10,20,50,80,100,120,150,180,200, 250, 300, 500, 1000, 2000, p['num_turns']-1] # turns to save particles object
p['turns2plot'] = [] # turns to plot phase space (while tracking)

# Beam intensity and emittance
p['n_part'] = int(1e4) # number of macroparticles
p['bunch_intensity'] = 40.0e10 # number of particles per bunch
p['macrosize'] = p['bunch_intensity']/p['n_part'] # number of charges per macroparticle
p['particle_distribution'] = 'real' # 'simulated' or 'real
if p['particle_distribution'] == 'simulated':
    p['nemitt_x'] = 0.6e-6 # normalized emittance in x
    p['nemitt_y'] = 0.6e-6 # normalized emittance in y
    p['sigma_z'] = (400/4)*0.525*0.3 # bunch length in m
    p['longitudinal_shape'] = 'gaussian' # 'parabolic' or 'coasting' or 'gaussian'
elif p['particle_distribution'] == 'real':
    p['nemitt_x'] = 0.4e-6 # guess, to be revised
    p['nemitt_y'] = 0.3e-6 # guess, to be revised
    p['sigma_z'] = 10 # guess, to be revised
    p['longitudinal_shape'] = np.nan # not used

# Tunes (at injection), chromaticity, imperfections
p['qx_ini'] = 4.40
p['qy_ini'] = 4.45
p['correct_chroma'] = False # if True, chromaticity is corrected using XNOH0
p['chroma_plane'] = 'y' # to correct horizontal ('x') or vertical ('y') chromaticity
p['include_field_errors'] = False # if True, field errors are included
p['field_errors'] = {
    'kbr1qno816l3': -6.15363e-4, # half-integer excitation (deltaI_816 = -2A)
    'kbr1qno412l3': 0,
}

# Setup acceleration
p['prepare_acceleration'] = 2 # 0: all RF OFF, 1: nominal PSB acceleration (double RF), 2: flat bottom (single RF at 8kV)
p['twiss_mode'] = '4d' # '4d' or '6d, used only if all RF are OFF
p['zeta0 = 17.5'] # if double RF, 6d-twiss at zeta0=0 fails because is unstable fixed point; this is a guess of the stable fixed point

# L4 parameters and number of injections
p['choppingFactor'] = 0.7
p['Linac4_current'] = 25e-3 # Amps
#p['num_injections'] = int(np.ceil(bunch_intensity/choppingFactor/Linac4_current/6.25e12))
p['num_injections'] = 1 # if > 1: a multi-turn injection is setup

# Injection chicane and correction
p['include_injection_chicane'] = 1 # if 1, 002A_include_injection_chicane.py is executed
p['on_chicane_k0 = 1'] # if 1, activates edge focusing of injection chicane
p['on_chicane_k2 = 1'] # if 1, activates eddy currents of injection chicane
p['include_injection_chicane_correction'] = 0 # if 1, 002B_include_injection_chicane_correction.py is executed
p['on_chicane_tune_corr = 1'] # if 1, activates tune correction of injection chicane
p['on_chicane_beta_corr = 1'] # if 1, activates beta correction of injection chicane

# Setup injection foil
p['install_injection_foil'] = False # if True, injection foil is installed
p['scatterchoice'] = 1 # 1: simple (no losses) 0: full (with losses) not working (for now)!

# Transverse painting
p['prepare_painting'] = 1 # if 1, 002C_prepare_painting.py is executed
p['on_painting_bump'] = 1 # if 1, activates painting bump
# Following convention: (t0,A0), (t1,A1), (t2,A2), (t3,A3), (1000,0)
# ISOLDE-like painting (check PSB_MD logbook entry 3806505 from 20/07/2023)
# https://logbook.cern.ch/elogbook-server/GET/showEventInLogbook/3806505
t2 = 148e-6 # num_injections*1e-6
t3 = 158e-6 # (num_injections+10)*1e-6
p['KSW_time_sec'] =         np.array([0.0e-6,  10e-6,      t2,     t3,  1000e-6])
p['KSW_bump_amplitude_m'] = np.array([-35e-3, -23e-3,  -12e-3, 9.2e-3,  0.0])

# Setup slicing and line cycling
p['slices'] = 3 # number of slices in thin lattice
# to have the starting point of the lattice at a different location, None otherwise
#p['element_to_cycle'] = None # line starts at the start of the 1st sector (NOT at the stripping foil)
p['element_to_cycle'] = 'bi1.tstr1l1' # stripping foil
#p['element_to_cycle'] = 'br1.bwsv11l1' # vertical LIU wire scanner

# Setup a linear tune ramp
p['prepare_tune_ramp'] = 0 # if 1, 004_prepare_tune_ramp.py is executed
p['on_tune_ramp'] = 1 # if 1, activates a linear tune ramp
p['qx_fin'] = 4.17 # final horizontal tune
p['qy_fin'] = 4.23 # final vertical tune

# Setup space charge calculation
p['install_space_charge'] = False # if True, space charge is installed
p['space_charge_mode'] = 'pic' # 'frozen' or 'pic' or 'quasi-frozen'
p['num_spacecharge_interactions'] = 160 # space charge interactions per turn
p['pic_solver'] = 'FFTSolver3D' # `FFTSolver2p5DAveraged` or `FFTSolver2p5D` or 'FFTSolver3D'

# Setup resources
p['GPU_FLAG'] = False # if True, GPU is used
if p['GPU_FLAG']:
    p['context'] = xo.ContextCupy()
else:
    p['context'] = xo.ContextCpu()

parameters = p