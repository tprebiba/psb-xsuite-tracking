import xtrack as xt
import pandas as pd
import xdeps as xd
import numpy as np
from simulation_parameters import parameters as p

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('003_prepare_acceleration.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

if p['prepare_acceleration']==0:
    print('prepare_acceleration = 0')
    print('No acceleration included.')

elif p['prepare_acceleration']==1:
    print('prepare_acceleration = 1')
    print('Double RF; with acceleration.')

    #########################################
    # Load PSB line in xsuite
    #########################################
    line = xt.Line.from_json('psb/psb_line_thick.json')
    line.build_tracker()

    #########################################
    # Load momentum and RF functions
    #########################################
    fname = 'time_tables/Ramp_and_RF_functions.dat'
    df = pd.read_csv(fname, sep='\t', skiprows=2,
        names=['t_s', 'E_kin_GeV', 'V1_MV', 'phi1_rad', 'V2_MV', 'phi2_rad'])
    print('Momentum and RF functions loaded from %s.'%(fname))

    #########################################
    # Set line energy program
    #########################################
    ep = xt.EnergyProgram(t_s=df['t_s'].values, kinetic_energy0=df['E_kin_GeV'].values*1e9)
    line.energy_program = ep
    print('Line energy program set.')

    ########################################
    # RF frequency functions
    ########################################
    t_s = df['t_s'].values
    freq_rev = line.energy_program.get_frev_at_t_s(t_s) # Revolution frequency at selected times
    #freq_rev

    # Create interpolating functions and variables
    line.functions['fun_freq_rev'] = xd.FunctionPieceWiseLinear(x=t_s, y=freq_rev)
    line.vars['freq_rev'] = line.functions['fun_freq_rev'](line.vars['t_turn_s'])

    # Create variables for first and second harmonics
    line.vars['freq_h1'] = line.vars['freq_rev']
    line.vars['freq_h2'] = 2 * line.vars['freq_rev']

    print('RF frequency functions built and attached to line.')

    ########################################
    # RF voltage and phase programs
    ########################################
    V1_MV = df.V1_MV.values # comes from the file
    V2_MV = df.V2_MV.values # comes from the file

    # Shift phases to have the beam centered around zero
    phi1_rad = df.phi1_rad.values - np.pi
    phi2_rad = df.phi2_rad.values - np.pi

    # Build interpolating functions and variables for voltages
    line.functions['fun_volt_mv_h1'] = xd.FunctionPieceWiseLinear(x=t_s, y=V1_MV)
    line.functions['fun_volt_mv_h2'] = xd.FunctionPieceWiseLinear(x=t_s, y=V2_MV)
    line.vars['volt_mv_h1'] = line.functions['fun_volt_mv_h1'](line.vars['t_turn_s'])
    line.vars['volt_mv_h2'] = line.functions['fun_volt_mv_h2'](line.vars['t_turn_s'])

    # Build interpolating functions and variables for voltages
    line.functions['fun_phi_rad_h1'] = xd.FunctionPieceWiseLinear(x=t_s, y=phi1_rad)
    line.functions['fun_phi_rad_h2'] = xd.FunctionPieceWiseLinear(x=t_s, y=phi2_rad)
    line.vars['phi_rad_h1'] = line.functions['fun_phi_rad_h1'](line.vars['t_turn_s'])
    line.vars['phi_rad_h2'] = line.functions['fun_phi_rad_h2'](line.vars['t_turn_s'])

    print('RF voltage and phase programs built and attached to line.')

    ########################################
    # Attach variables to cavities
    ########################################
    line.element_refs['br1.acwf5l1.1'].voltage = line.vars['volt_mv_h1'] * 1e6
    line.element_refs['br1.acwf5l1.2'].voltage = line.vars['volt_mv_h2'] * 1e6
    line.element_refs['br1.acwf5l1.1'].lag = line.vars['phi_rad_h1'] * 360 / 2 / np.pi
    line.element_refs['br1.acwf5l1.2'].lag = line.vars['phi_rad_h2'] * 360 / 2 / np.pi
    line.element_refs['br1.acwf5l1.1'].frequency = line.vars['freq_h1']
    line.element_refs['br1.acwf5l1.2'].frequency = line.vars['freq_h2']
    print('Voltage and phase functions attached to cavities.')

    #########################################
    # Save line to .json
    #########################################
    line.vars['t_turn_s']=0.0 # Reset time at the injection
    #line.twiss_default['method'] = '6d' # now with cavity
    #print('Twiss method set to 6d.')
    line.twiss()
    line.to_json('psb/psb_line_thick.json')
    print('Line saved to psb/psb_line_thick.json')

elif p['prepare_acceleration']==2: 
    print('prepare_acceleration = 2')
    print('Single RF; no acceleration included.')

    #########################################
    # Load PSB line in xsuite
    #########################################
    line = xt.Line.from_json('psb/psb_line_thick.json')
    line.build_tracker()

    #########################################
    # Add constant voltage to the dummy RF
    #########################################
    line.element_refs['br.c02'].voltage = 0.008*1e6
    print('Constant voltage = 8 kV added to the dummy RF.')

    #########################################
    # Save line to .json
    #########################################
    line.vars['t_turn_s']=0.0 # Reset time at the injection
    line.twiss_default['method'] = '6d' # now with cavity
    print('Twiss method set to 6d.')
    line.twiss()
    line.to_json('psb/psb_line_thick.json')
    print('Line saved to psb/psb_line_thick.json')