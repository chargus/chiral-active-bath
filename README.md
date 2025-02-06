# LAMMPS simulations of a chiral active bath
![Chiral active bath](cover_figure.png)
Simulations of passive objects in a chiral active bath

Simulation and analysis code.

To reproduce the simulations of passive objects in a bath of chiral ABPs, as described in [C. Hargus, F. Ghimenti, J. Tailleur and F. van Wijland arXiv:2412.20689 (2025)], you should take the following steps:
    
1. Download LAMMPS from https://lammps.sandia.gov/
2. Add the custom source and header files contained in this repository (in ./src) to the LAMMPS src directory (e.g. /path/to/lammps/src)
3. Compile LAMMPS following instructions at https://lammps.sandia.gov/
4. Invoke the lmp binary on a script from ./scripts. For example, `mpirun -np 24 lmp -in disk_harmonic.in`

