import xtrack as xt
from simulation_parameters import parameters as p

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')
print('006_lattice_imperfections.py')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~**~*~*~**~*~*~*')

#########################################
# Load PSB line in xsuite
#########################################
line = xt.Line.from_json('psb/psb_line_thin.json') # to get the lengths of the injection bumpers

#########################################
# Correct vertical chromaticity
#########################################
if p['correct_chroma']:
    print('Correcting chromaticity in %s plane.'%(p['chroma_plane']))
    line.match(vary=[xt.Vary('kbr1xnoh0', step=1e-8)],
               targets = [xt.Target('dq%s'%(p['chroma_plane']), 0.0, tol=1e-5)])
    print('Chromaticity corrected.')
else:
    print('Natural chromaticity.')

#########################################
# Add field errors
#########################################
if p['include_field_errors']:
    print('Adding field errors.')
    for key in p['field_errors']:
        line.vars[key] = p['field_errors'][key]
        print('Added %s = %s to the lattice.'%(key, p['field_errors'][key]))
else:
    print('No field errors.')

#########################################
# Save line to .json
#########################################
line.to_json('psb/psb_line_thin.json')
print('Line saved to psb/psb_line_thin.json')