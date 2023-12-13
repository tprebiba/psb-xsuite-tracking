import xtrack as xt
import xdeps as xd
import numpy as np
from simulation_parameters import parameters as p

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('002B_include_injection_chicane_correction.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

if p['include_injection_chicane_correction']==0:
    print('include_injection_chicane_correction = 0')
    print('Injection correction will not be included.')

elif p['include_injection_chicane_correction']==1:

    #########################################
    # Load PSB line in xsuite
    #########################################
    line = xt.Line.from_json('psb/psb_line_thick.json')
    line.build_tracker()


    #########################################
    # Insert marker at location where
    # we want to mimize beta-beating
    #########################################
    line.discard_tracker() # We need to discard the tracker to edit the line
    line.insert_element(element=xt.Marker(), name='mker_match', at_s=79.874)
    line.build_tracker()


    #########################################
    # Store twiss with chicane off
    #########################################
    line.vars['on_chicane_k0'] = 0
    line.vars['on_chicane_k2'] = 0
    tw0 = line.twiss()
    line.vars['on_chicane_k0'] = 1
    line.vars['on_chicane_k2'] = 1


    #########################################
    # Compute correction using xsuite match
    #########################################
    print('Computing correction using xsuite match...')

    t_correct = np.linspace(0, 5.5e-3, 30) # Times at which corrections are computed

    kbrqf_corr_list = []
    kbrqd_corr_list = []
    kbrqd3corr_list = []
    kbrqd14corr_list = []

    for ii, tt in enumerate(t_correct):
        print(f'Correct tune and beta beat at t = {tt * 1e3:.2f} ms', end='\r', flush=True)
        
        line.vars['t_turn_s'] = tt # Set test time

        xd.general._print.suppress = True # to avoid too much output
        line.match(
            targets=[
                xt.Target('qx', value=tw0.qx, tol=1e-5, scale=1),
                xt.Target('qy', value=tw0.qy, tol=1e-5, scale=1),
                xt.Target('bety', at='mker_match',
                        value=tw0['bety', 'mker_match'], tol=1e-4, scale=100),
                xt.Target('alfy', at='mker_match',
                        value=tw0['alfy', 'mker_match'], tol=1e-4, scale=100)
            ],
            vary=[
                xt.Vary('kbrqfcorr', step=1e-4),
                xt.Vary('kbrqdcorr', step=1e-4),
                xt.Vary('kbrqd3corr', step=1e-4),
                xt.Vary('kbrqd14corr', step=1e-4),
            ],
        )
        xd.general._print.suppress = False

        # Store found values
        kbrqf_corr_list.append(line.vars['kbrqfcorr']._value)
        kbrqd_corr_list.append(line.vars['kbrqdcorr']._value)
        kbrqd3corr_list.append(line.vars['kbrqd3corr']._value)
        kbrqd14corr_list.append(line.vars['kbrqd14corr']._value)
    line.vars['t_turn_s']=0.0 # Reset time at the injection
    print('Correction computed using xsuite match.')


    #########################################
    # Build function with the computed 
    # corrections
    #########################################
    line.functions['fun_kqf_corr'] = xd.FunctionPieceWiseLinear(x=t_correct, y=kbrqf_corr_list)
    line.functions['fun_kqd_corr'] = xd.FunctionPieceWiseLinear(x=t_correct, y=kbrqd_corr_list)
    line.functions['fun_qd3_corr'] = xd.FunctionPieceWiseLinear(x=t_correct, y=kbrqd3corr_list)
    line.functions['fun_qd14_corr'] = xd.FunctionPieceWiseLinear(x=t_correct, y=kbrqd14corr_list)
    print('Beta-beating correction functions built to control quad strengths.')


    #########################################
    # Use functions to control quad strengths
    #########################################
    line.vars['on_chicane_tune_corr'] = 1
    line.vars['kbrqfcorr'] = (line.vars['on_chicane_tune_corr']
                                * line.functions.fun_kqf_corr(line.vars['t_turn_s']))
    line.vars['kbrqdcorr'] = (line.vars['on_chicane_tune_corr']
                                * line.functions.fun_kqd_corr(line.vars['t_turn_s']))

    line.vars['on_chicane_beta_corr'] = 1
    line.vars['kbrqd3corr'] = (line.vars['on_chicane_beta_corr']
                            * line.functions.fun_qd3_corr(line.vars['t_turn_s']))
    line.vars['kbrqd14corr'] = (line.vars['on_chicane_beta_corr']
                            * line.functions.fun_qd14_corr(line.vars['t_turn_s']))
    print('Variables "on_chicane_tune_corr" and "on_chicane_beta_corr" control the scaling of the tune & beta-beating correction.')
    print('Beta-beating correction function assigned to quadrupoles.')


    #########################################
    # Save line to .json
    #########################################
    line.to_json('psb/psb_line_thick.json')
    print('Line saved to psb/psb_line_thick.json')