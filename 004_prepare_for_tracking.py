import xtrack as xt
from xtrack.slicing import Strategy, Teapot
from simulation_parameters import parameters as p

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('004_prepare_for_tracking.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

if p['prepare_tune_ramp']==0:
    print('prepare_tune_ramp = 0')
    print('Tunes will remain constant.')

elif p['prepare_tune_ramp']==1:

    #########################################
    # Load PSB line in xsuite
    #########################################
    line = xt.Line.from_json('psb/psb_line_thick.json') # to get the lengths of the injection bumpers

    #########################################
    # Deactivate chicane and correction
    #########################################
    line.vars['on_chicane_k0'] = 0
    line.vars['on_chicane_k2'] = 0
    line.vars['on_chicane_tune_corr'] = 0
    line.vars['on_chicane_beta_corr'] = 0
    print('Chicane and correction deactivated.')

    #########################################
    # Slice thick line with xsuite 
    #########################################
    slices = 3
    line.slice_thick_elements(
        slicing_strategies=[
            Strategy(slicing=Teapot(slices)), # default is 1
            Strategy(slicing=Teapot(slices), element_type=xt.Bend),
            Strategy(slicing=Teapot(slices), element_type=xt.Quadrupole),
            Strategy(slicing=Teapot(slices), name=r'bi.*bsw.*'),
        ]
    )
    print('Line sliced into %s slices with xsuite.'%(slices))

    #########################################
    # Rematch tunes in thin lattice and twiss
    #########################################
    print('Matching with xsuite.')
    qx_target = p['qx_ini']
    qy_target = p['qy_ini']
    line.match(
                vary=[
                        xt.Vary('kbrqf', step=1e-8),
                        xt.Vary('kbrqd', step=1e-8),
                ],
                targets = [
                            xt.Target('qx', qx_target, tol=1e-5),
                            xt.Target('qy', qy_target, tol=1e-5)
                ]
    )
    tw = line.twiss()
    print('Twiss computed.')
    print('Working point of thin lattice: (Qx, Qy) = (%s, %s)'%(tw.qx, tw.qy))

    #########################################
    # Reactivate chicane and correction
    #########################################
    line.vars['on_chicane_k0'] = p['on_chicane_k0']
    line.vars['on_chicane_k2'] = p['on_chicane_k2']
    line.vars['on_chicane_tune_corr'] = p['on_chicane_tune_corr']
    line.vars['on_chicane_beta_corr'] = p['on_chicane_beta_corr']
    print('Chicane and correction reactivated.')
    tw = line.twiss()
    tw.to_pandas().to_pickle('psb/psb_twiss_thin.pkl')
    print('Twiss computed and saved to psb/psb_twiss_thin.pkl.')
    print('Working point of thin lattice: (Qx, Qy) = (%s, %s)'%(tw.qx, tw.qy))

    #########################################
    # Start lattice at desired location
    #########################################
    line.discard_tracker() # We need to discard the tracker to edit the line
    line.cycle(name_first_element = 'br1.bwsv11l1', inplace=True) # vertical (LIU) wire scanner
    print('Changed line starting point to br1.bwsv11l1 (vertical wire scanner).')

    #########################################
    # Save line to .json
    #########################################
    line.to_json('psb/psb_line_thin.json')
    print('Line saved to psb/psb_line_thin.json')