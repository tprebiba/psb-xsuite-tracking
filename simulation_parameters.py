#######################################
# Simulation parameters
#######################################
n_part = int(5e5) #int(5e5) # number of macroparticles: 2e5 were used in pyorbit
#n_part = int(4e5)

num_turns = 100000 #100000 # number of turns to track

bunch_intensity = 40.0e10 # number of particles per bunch

nemitt_x = 1.05e-6 # normalized emittance in x
nemitt_y = 0.62e-6 # normalized emittance in y
sigma_z = (400/4)*0.525*0.3 # bunch length in m ?

num_spacecharge_interactions = 160 # space charge interactions per turn
tol_spacecharge_position = 1e-2 # minimum/maximum space between sc elements

qx_ini = 4.15
qy_ini = 4.45

# Half-integer excitation (deltaI_816 = -2A)
kbr1qno816l3 = -6.15363e-4#*10
kbr1qno412l3 = 0


#######################################
# Flags to control simulation flow
#######################################
include_injection_chicane = 0 # if 1, 002A_include_injection_chicane.py is executed
on_chicane_k0 = 1 # if 1, activates edge focusing of injection chicane
on_chicane_k2 = 1 # if 1, activates eddy currents of injection chicane

include_injection_chicane_correction = 0 # if 1, 002B_include_injection_chicane_correction.py is executed
on_chicane_tune_corr = 1 # if 1, activates tune correction of injection chicane
on_chicane_beta_corr = 1 # if 1, activates beta correction of injection chicane

prepare_acceleration = 0 # if 1, 003_prepare_acceleration.py is executed

prepare_tune_ramp = 0 # if 1, 004_prepare_tune_ramp.py is executed
on_tune_ramp = 1 # if 1, activates tune ramp
qx_fin = 4.17 # final horizontal tune
qy_fin = 4.23 # final vertical tune


#######################################
# Store parameters in dictionary
#######################################
parameters = {
    'n_part': n_part,
    'num_turns': num_turns,
    'bunch_intensity': bunch_intensity,
    'nemitt_x': nemitt_x,
    'nemitt_y': nemitt_y,
    'sigma_z': sigma_z,
    'num_spacecharge_interactions': num_spacecharge_interactions,
    'tol_spacecharge_position': tol_spacecharge_position,
    'qx_ini': qx_ini,
    'qy_ini': qy_ini,
    'kbr1qno816l3': kbr1qno816l3,
    'kbr1qno412l3': kbr1qno412l3,
    'include_injection_chicane': include_injection_chicane,
    'on_chicane_k0': on_chicane_k0,
    'on_chicane_k2': on_chicane_k2,
    'include_injection_chicane_correction': include_injection_chicane_correction,
    'on_chicane_tune_corr': on_chicane_tune_corr,
    'on_chicane_beta_corr': on_chicane_beta_corr,
    'prepare_acceleration': prepare_acceleration,
    'prepare_tune_ramp': prepare_tune_ramp,
    'on_tune_ramp': on_tune_ramp,
    'qx_fin': qx_fin,
    'qy_fin': qy_fin,
}
