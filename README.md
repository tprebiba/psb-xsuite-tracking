# Realistic Proton Synchrotron Booster (PSB) Beam Dynamics Simulation with Xsuite

This repository provides a full simulation pipeline for realistic modeling of beam dynamics in CERN’s Proton Synchrotron Booster (PSB) using [Xsuite](https://xsuite.readthedocs.io/en/latest/).  

It includes tools for:
- Lattice generation
- Multi-turn injection and transverse painting modeling  
- Foil scattering
- Acceleration
- Time-dependent settings (e.g. tune ramps)
- Machine imperfections  
- Space charge effects  
- Tracking and output analysis  
- Scripts for submittion to [HTCondor](https://abpcomputing.web.cern.ch/computing_resources/cernbatch/) (CPU or GPU)

**Contact for corrections & suggestions**: [tirsi.prebibaj@cern.ch](mailto:tirsi.prebibaj@cern.ch)  
With inputs from: F. Asvesta, H. Bartosik, G. Iadarola, K. Paraschou

---

## Structure

The simulation is organized in two main parts:

### Part I: Preparing the Simulation

- Generate the desired lattice and particle distribution (lattice setup, beam transverse and longitudinal characteristics, ...).
- All settings are controlled via `simulation_parameters.py`.
- Execution:
    ```
    . 000_prepare_simulation.sh 
    ```
- Will generate the lattice and machine settings in `psb/` and the initial particle distribution in `input/`


### Part II: Tracking

- Perform beam tracking (configured also via `simulation_parameters.py`).
- Local execution:
    ```
    python -m runPSB.py 
    ```
    or execution in HTCondor with GPU
    ```
    condor_submit htcondor_submission_gpu.sub
    ```
- All the outputs (turn-by-turn beam data, beam profiles, ...) are saved in `output/`.