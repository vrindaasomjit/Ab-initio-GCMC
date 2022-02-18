#!/bin/bash
#SBATCH -J gcmc             # Job Name
#SBATCH -o gcmc.o%j             # Output file name (%j expands to jobID)
#SBATCH -e gcmc.e%j             # Error file name (%j expands to jobID)
#SBATCH -N 4                    # Number of nodes
#SBATCH -n 182
#SBATCH -p skx-dev               # Queue name (skx-normal, development, etc.)
#SBATCH -t 2:00:00             # Run time (hh:mm:ss)

module reset
module load intel/18.0.2
module load impi/18.0.2
module load python2/2.7.15
module load vasp/5.4.4
module unload mistral/2.13.4
ml

export LAMMPS_DIR=/home1/04770/tg840694/help_TACC_lammps/stable_3Mar2020_clean/
export PATH=${PATH}:${LAMMPS_DIR}/bin
which lmp_stampede

export PYTHONPATH=${PYTHONPATH}:${LAMMPS_DIR}/lib/message/cslib/src
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LAMMPS_DIR}/lib/message/cslib/src

#ibrun -n 1 -o 0 lmp_stampede -restart2data gcmc.restart.1353 remap data.intf  
ibrun -n 1 -o 0 lmp_stampede -v mode file < in.client.intf &
python vasp_wrap_gcmc.py file POSCARintf &
wait

