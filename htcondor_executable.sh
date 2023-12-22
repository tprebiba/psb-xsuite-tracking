#!/bin/bash
echo "uname -r:" `uname -r`

#export MYPYTHON=/afs/cern.ch/user/t/tprebiba/python_environments/base-2022/miniconda
source /usr/local/xsuite/miniforge3/bin/activate xsuite
pwd

current_directory="/home/tprebiba/eos/Fellowship/03_Xsuite/psb-xsuite-tracking"

cp -r $current_directory/psb .
cp -r $current_directory/input .
cp -r $current_directory/output .
cp -r $current_directory/lib .
cp -r $current_directory/time_tables .
cp $current_directory/*.py .
cp $current_directory/*.sh .

#source $MYPYTHON/bin/activate ""
#export PATH=$MYPYTHON/bin:$PATH
echo "which python:" `which python`

date
#DIR=/afs/cern.ch/work/t/tprebiba/half-integer_xsuite/000_PSB_half-integer_dynamic_crossing_PIC_reference
#python $DIR/runPSB.py
python runPSB.py
date
