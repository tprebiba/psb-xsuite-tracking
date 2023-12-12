#!/bin/bash
echo "uname -r:" `uname -r`

#export MYPYTHON=/afs/cern.ch/user/t/tprebiba/python_environments/base-2022/miniconda
source /usr/local/xsuite/miniforge3/bin/activate xsuite
pwd

cp -r /afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference/psb .
cp -r /afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference/tables .
cp -r /afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference/input .
cp -r /afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference/output .
cp -r /afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference/lib .
cp /afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference/*.py .
cp /afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference/*.sh .

#source $MYPYTHON/bin/activate ""
#export PATH=$MYPYTHON/bin:$PATH
echo "which python:" `which python`

date
#DIR=/afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference
#python $DIR/runPSB.py
python runPSB.py
date
