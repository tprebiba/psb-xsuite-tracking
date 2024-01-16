**Realistic PSB simulation in [Xsuite](https://xsuite.readthedocs.io/en/latest/).**

For corrections & suggestions contact: [tirsi.prebibaj@cern.ch](mailto:tirsi.prebibaj@cern.ch)
With inputs from F. Asvesta, H. Bartosik, G. Iadarola, K. Paraschou. 

Simulation is split in two parts:

1. Part I: preparing the simulation
   - Generate the desired lattice and particle distribution (set tune, injection bump, foil scattering, acceleration, imperfections, time-functions, emittances, bunch shape, …).​
   - Inputs are given by ```simulation_parameters.py```.
   - Execution:
       ```
        $ . 000_prepare_simulation.sh
       ```
   - Outputs stored in ```psb/``` and ```input/```.
   - Runs locally on CPU

3. Part II: tracking
   - Track and save beam state (space charge algorithms, multi-turn injection, …).​
   - Inputs are given by ```simulation_parameters.py```.
   - Execution:
       ```
        $ python -m runPSB.py
       ```
   - Outputs stored in ```output/```.
   - Runs locally or on HTCondor, in CPU or GPU.
   - Run on HTCondor using GPUs with:
        ```
        $ condor_submit htcondor_submission_gpu.sub
       ```