import xtrack as xt
import xdeps as xd
import pandas as pd
from simulation_parameters import parameters as p

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('002C_prepare_painting.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

if p['prepare_painting']==0:
    print('prepare_painting = 0')
    print('No painting applied.')

elif p['prepare_painting']==1:

    #########################################
    # Load PSB line in xsuite
    #########################################
    line = xt.Line.from_json('psb/psb_line_thick.json')
    line.build_tracker()

    #########################################
    # Load KSW bump file
    #########################################
    fname = 'time_tables/KSW.tfs'
    df = pd.read_csv(fname, sep='\s+', skiprows=8,
        #names=['X_BUMP', 'KBIKSW1L4', 'KBIKSW2L1', 'KBIKSW16L1', 'KBIKSW16L4'])
        names=['X_BUMP', 'kbi1ksw1l4', 'kbi1ksw2l1', 'kbi1ksw16l1', 'kbi1ksw16l4']) # different naming convention
    print('KSW bump file loaded.')

    #########################################
    # A few checks
    #########################################
    # tw = line.twiss()
    # co_x_at_foil = tw['x', 'bi1.tstr1l1_entry']
    # co_y_at_foil = tw['y', 'bi1.tstr1l1_entry']
    # print('Closed orbit at foil: x = %s m, y = %s m.'%(co_x_at_foil, co_y_at_foil))
    # line.vars['kbi1ksw1l4'] = df['kbi1ksw1l4'].values[0]*.37
    # line.vars['kbi1ksw2l1'] = df['kbi1ksw2l1'].values[0]*.37
    # line.vars['kbi1ksw16l1'] = df['kbi1ksw16l1'].values[0]*.37
    # line.vars['kbi1ksw16l4'] = df['kbi1ksw16l4'].values[0]*.37
    # tw = line.twiss()
    # co_x_at_foil = tw['x', 'bi1.tstr1l1_entry']
    # co_y_at_foil = tw['y', 'bi1.tstr1l1_entry']
    # print('Closed orbit at foil: x = %s m, y = %s m.'%(co_x_at_foil, co_y_at_foil)) # should be x=+35mm

    #########################################
    # Prepare knobs and add time dependence
    #########################################
    # line.element_dict['bi1.ksw1l4'].length is zero although it should be .37 m
    # hardcoded for now, to be rechecked
    ksw_length = 0.37 # m
    # strengths need a -1 factor to produce a -0.035 m bump, to be rechecked
    line.vars['on_painting_bump'] = p['on_painting_bump']
    kbi1ksw1l4 = (-1)*p['KSW_bump_amplitude_m']/df['X_BUMP'].values[0]*df['kbi1ksw1l4'].values[0]*ksw_length
    kbi1ksw2l1 = (-1)*p['KSW_bump_amplitude_m']/df['X_BUMP'].values[0]*df['kbi1ksw2l1'].values[0]*ksw_length
    kbi1ksw16l1 = (-1)*p['KSW_bump_amplitude_m']/df['X_BUMP'].values[0]*df['kbi1ksw16l1'].values[0]*ksw_length
    kbi1ksw16l4 = (-1)*p['KSW_bump_amplitude_m']/df['X_BUMP'].values[0]*df['kbi1ksw16l4'].values[0]*ksw_length
    line.functions['fun_kbiksw1l4'] = xd.FunctionPieceWiseLinear(x=p['KSW_time_sec'], y=kbi1ksw1l4)
    line.functions['fun_kbiksw2l1'] = xd.FunctionPieceWiseLinear(x=p['KSW_time_sec'], y=kbi1ksw2l1)
    line.functions['fun_kbiksw16l1'] = xd.FunctionPieceWiseLinear(x=p['KSW_time_sec'], y=kbi1ksw16l1)
    line.functions['fun_kbiksw16l4'] = xd.FunctionPieceWiseLinear(x=p['KSW_time_sec'], y=kbi1ksw16l4)
    line.vars['kbi1ksw1l4'] = line.functions.fun_kbiksw1l4(line.vars['t_turn_s'])*line.vars['on_painting_bump']
    line.vars['kbi1ksw2l1'] = line.functions.fun_kbiksw2l1(line.vars['t_turn_s'])*line.vars['on_painting_bump']
    line.vars['kbi1ksw16l1'] = line.functions.fun_kbiksw16l1(line.vars['t_turn_s'])*line.vars['on_painting_bump']
    line.vars['kbi1ksw16l4'] = line.functions.fun_kbiksw16l4(line.vars['t_turn_s'])*line.vars['on_painting_bump']
    print('Knobs built to control KSW bump.')
    print('Variable "on_painting_bump" added to control KSW bump.')
    print('Time dependence added to kbi1ksw1l4, kbi1ksw2l1, kbi1ksw16l1 and kbi1ksw16l4.')

    #########################################
    # Save line to .json
    #########################################
    line.to_json('psb/psb_line_thick.json')
    print('Line saved to psb/psb_line_thick.json')