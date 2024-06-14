import xtrack as xt
from xtrack.slicing import Strategy, Teapot
from simulation_parameters import parameters as p

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('004_prepare_for_tracking.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

#########################################
# Load PSB line in xsuite
#########################################
line = xt.Line.from_json('psb/psb_line_thick.json') # to get the lengths of the injection bumpers

#########################################
# Deactivate chicane and correction
# Deactivate painting bump
#########################################
line.vars['on_chicane_k0'] = 0
line.vars['on_chicane_k2'] = 0
line.vars['on_chicane_tune_corr'] = 0
line.vars['on_chicane_beta_corr'] = 0
line.vars['on_painting_bump'] = 0
print('Chicane and correction deactivated.')
print('Painting bump deactivated.')

#########################################
# Slice thick line with xsuite 
#########################################
slices = p['slices']
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
# Reactivate painting bump
#########################################
line.vars['on_chicane_k0'] = p['on_chicane_k0']
line.vars['on_chicane_k2'] = p['on_chicane_k2']
line.vars['on_chicane_tune_corr'] = p['on_chicane_tune_corr']
line.vars['on_chicane_beta_corr'] = p['on_chicane_beta_corr']
line.vars['on_painting_bump'] = p['on_painting_bump']
if ((p['include_injection_chicane']>0) or (p['include_injection_chicane_correction']>0)):
    print('Chicane and correction reactivated.')
if p['prepare_painting']>0:
    print('Painting bump reactivated.')
tw = line.twiss()
tw.to_pandas().to_pickle('psb/psb_twiss_thin.pkl')
print('Twiss computed and saved to psb/psb_twiss_thin.pkl.')
print('Working point of thin lattice: (Qx, Qy) = (%s, %s)'%(tw.qx, tw.qy))

#########################################
# Save line to .json
#########################################
line.to_json('psb/psb_line_thin.json')
print('Line saved to psb/psb_line_thin.json')