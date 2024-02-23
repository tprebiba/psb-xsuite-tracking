#!/bin/bash

# This script prepares the simulation directory for the simulation

# Pull psb lattice from acc-models-psb repository on gitlab.cern.ch
pull_acc_models_psb_files=false
if [ "$pull_acc_models_psb_files" = true ]; then
    echo "Pulling psb lattice from acc-models-psb repository on gitlab.cern.ch"
    cd psb
    echo "Current directory: $(pwd)"
    source pull_acc-models-psb_files.sh
    echo "Pulled."
    cd ..
    echo "Current directory: $(pwd)"
fi

# Configure madx lattice for xsuite (matching, cycling, etc.)
python 001_get_PSB_line.py

# Configure injection chicane
python 002A_include_injection_chicane.py

# Configure injection chicane correction functions
python 002B_include_injection_chicane_correction.py

# Configure transverse painting
python 002C_prepare_painting.py

# prepare acceleration
python 003_prepare_acceleration.py

# Configure xsuite line for tracking (slicing, re-matching, etc.)
python 004_prepare_for_tracking.py

# Configure tune ramp
python 005_prepare_tune_ramp.py

# Include lattice imperfections (chroma, field errors, etc.)
python 006_lattice_imperfections.py

# Particle distribution (multi-turn injection)
python 007_generate_particle_distribution.py

# Configure paths on htcondor_executable for submission to htcondor
echo "Configuring paths on htcondor_executable for submission to htcondor"
current_directory=$(pwd)
htcondor_executable_file="htcondor_executable.sh"
if [ -e "$htcondor_executable_file" ]; then
    sed -i "8s|.*|current_directory=$current_directory|" "$htcondor_executable_file"
    echo "Current directory $current_directory written to $htcondor_executable_file."
else
    echo "Error: $htcondor_executable_file does not exist."
fi