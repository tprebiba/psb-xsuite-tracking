import xtrack as xt
import json
import xdeps as xd
import numpy as np
from simulation_parameters import parameters as p

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('005_prepare_tune_ramp.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

if p['prepare_tune_ramp']==0:
    print('prepare_tune_ramp = 0')
    print('Tunes will remain constant.')

elif p['prepare_tune_ramp']==1:

      #########################################
      # Load PSB line in xsuite
      #########################################
      line = xt.Line.from_json('psb/psb_line_thin.json')
      line.build_tracker()


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
      # Generate tune ramp
      #########################################
      print('Matching at (Qx, Qy) = (%s, %s)'%(p['qx_fin'], p['qy_fin']))
      qx_target_fin = p['qx_fin']
      qy_target_fin = p['qy_fin']
      line.match(
                  vary=[
                        xt.Vary('kbrqf', step=1e-8),
                        xt.Vary('kbrqd', step=1e-8),
                  ],
                  targets = [
                              xt.Target('qx', qx_target_fin, tol=1e-5),
                              xt.Target('qy', qy_target_fin, tol=1e-5)
                  ]
      )
      kf_fin = line.vars['kbrqf']._value
      kd_fin = line.vars['kbrqd']._value
      print('Converged to strengths (kf, kd) = (%s, %s)'%(kf_fin, kd_fin))

      print('Matching at (Qx, Qy) = (%s, %s)'%(p['qx_ini'], p['qy_ini']))
      qx_target_ini = p['qx_ini']
      qy_target_ini = p['qy_ini']
      line.match(
                  vary=[
                        xt.Vary('kbrqf', step=1e-8),
                        xt.Vary('kbrqd', step=1e-8),
                  ],
                  targets = [
                              xt.Target('qx', qx_target_ini, tol=1e-5),
                              xt.Target('qy', qy_target_ini, tol=1e-5)
                  ]
      )
      kf_ini = line.vars['kbrqf']._value
      kd_ini = line.vars['kbrqd']._value
      print('Converged to strengths (kf, kd) = (%s, %s)'%(kf_ini, kd_ini))

      d = {
            'qx_ini': p['qx_ini'],
            'qy_ini': p['qy_ini'],
            'kf_ini': kf_ini,
            'kd_ini': kd_ini,
            'qx_fin': p['qx_fin'],
            'qy_fin': p['qy_fin'],
            'kf_fin': kf_fin,
            'kd_fin': kd_fin,
            'num_turns': p['num_turns'],
      }
      with open('time_tables/tunes.json','w') as fid:
            json.dump(d, fid, indent=2)
      print('Dictionary: ', d)
      print('Wrote time_tables/tunes.json')


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
      print('Working point of thin lattice: (Qx, Qy) = (%s, %s)'%(tw.qx, tw.qy))


      #########################################
      # Generating xsuite functions and
      # assigning to main quads
      #########################################
      line.functions['kbrqf_func'] = xd.FunctionPieceWiseLinear(x=np.array([0, p['num_turns']/1e6]), y = np.array([kf_ini, kf_fin]))
      line.functions['kbrqd_func'] = xd.FunctionPieceWiseLinear(x=np.array([0, p['num_turns']/1e6]), y = np.array([kd_ini, kd_fin]))

      line.vars['on_tune_ramp'] = p['on_tune_ramp'] # to easily switch off the tune ramp
      line.vars['kbrqf'] = line.functions.kbrqf_func(line.vars['t_turn_s']) * line.vars['on_tune_ramp']
      line.vars['kbrqd'] = line.functions.kbrqd_func(line.vars['t_turn_s']) * line.vars['on_tune_ramp'] 


      #########################################
      # Save line to .json
      #########################################
      line.to_json('psb/psb_line_thin.json')
      print('Line saved to psb/psb_line_thin.json')