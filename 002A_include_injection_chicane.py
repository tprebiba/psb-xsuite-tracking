import xtrack as xt
import pandas as pd
import xdeps as xd
from simulation_parameters import parameters as p

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('002_include_injection_chicane.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

if p['include_injection_chicane']==0:
    print('include_injection_chicane = 0')
    print('Injection chicane collapse will not be included.')

elif p['include_injection_chicane']==1:

    #########################################
    # Load PSB line in xsuite
    #########################################
    line = xt.Line.from_json('psb/psb_line_thick.json')
    line.build_tracker()


    #########################################
    # A few checks on the imported model
    #########################################
    #line.vars['k0bi1bsw1l11']._info() # Check that the knob controls k0 and the edges
    #line.element_refs['bi1.bsw1l1.1'].h._info() # Check no reference system curvature in bumpers


    #########################################
    # Build knobs to model k0 and k2 
    # (edge focusing and eddy currents)
    #########################################
    line.vars['bsw_k0l'] = 0 # knob controlling the four magnets
    line.vars['k0bi1bsw1l11'] =  line.vars['bsw_k0l'] / line['bi1.bsw1l1.1'].length
    line.vars['k0bi1bsw1l12'] = -line.vars['bsw_k0l'] / line['bi1.bsw1l1.2'].length
    line.vars['k0bi1bsw1l13'] = -line.vars['bsw_k0l'] / line['bi1.bsw1l1.3'].length
    line.vars['k0bi1bsw1l14'] =  line.vars['bsw_k0l'] / line['bi1.bsw1l1.4'].length
    line.vars['bsw_k2l'] = 0 # knob controlling the four magnets
    line.element_refs['bi1.bsw1l1.1'].knl[2] = line.vars['bsw_k2l']
    line.element_refs['bi1.bsw1l1.2'].knl[2] = -line.vars['bsw_k2l']
    line.element_refs['bi1.bsw1l1.3'].knl[2] = -line.vars['bsw_k2l']
    line.element_refs['bi1.bsw1l1.4'].knl[2] = line.vars['bsw_k2l']
    print('Knobs built to control k0 and k2.')


    #########################################
    # Add time dependence
    #########################################
    df = pd.read_csv('time_tables/chicane_collapse.csv', delimiter=',', skipinitialspace=True)
    #df.head()
    line.functions['fun_bsw_k0l'] = xd.FunctionPieceWiseLinear(x=df['time'].values, y=df['bsw_k0l'].values)
    line.functions['fun_bsw_k2l'] = xd.FunctionPieceWiseLinear(x=df['time'].values, y=df['bsw_k2l'].values)
    line.vars['on_chicane_k0'] = p['on_chicane_k0'] # to easily switch off the effect of the edge focusing
    line.vars['on_chicane_k2'] = p['on_chicane_k2'] # to easily switch off the effect of the eddy currents
    line.vars['bsw_k0l'] = line.functions.fun_bsw_k0l(line.vars['t_turn_s']) * line.vars['on_chicane_k0']
    line.vars['bsw_k2l'] = line.functions.fun_bsw_k2l(line.vars['t_turn_s']) * line.vars['on_chicane_k2'] 
    print('Variables "on_chicane_k0" and "on_chicane_k2" added to switch off edge focusing and eddy currents.')
    print('Time dependence added to k0 and k2.')


    #########################################
    # Save line to .json
    #########################################
    line.vars['t_turn_s']=0.0 # Reset time at the injection
    line.to_json('psb/psb_line_thick.json')
    print('Line saved to psb/psb_line_thick.json')