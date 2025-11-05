#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Stampede2 KNL nodes
#
#   *** OpenMP Job on Normal Queue ***

#SBATCH -J mergertree           # Job name
#SBATCH -o myjob.oNUM       # Name of stdout output file
#SBATCH -e myjob.eNUM     # Name of stderr error file
#SBATCH -p skx          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 1              # Total # of mpi tasks
#SBATCH -t 4:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE    # Send email at begin and end of job
#SBATCH -A TG-AST140041

# Other commands must follow all #SBATCH directives...

#module load gcc
module load gcc
module load python
module list
pwd
date

#Set thread count (default value is 1)...

#export OMP_NUM_THREADS=48


# Launch OpenMP code...
rm output.txt
python -u find_leo_candidates.py START STOP >> outputNUM.txt

date
