#!/bin/bash
#SBATCH -J mergertree           # Job name
#SBATCH -o myjob.oNUM     # Name of stdout output file
#SBATCH -e myjob.eNUM      # Name of stderr error file
#SBATCH -p skx          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 1              # Total # of mpi tasks
#SBATCH -t 2:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE    # Send email at begin and end of job
#SBATCH -A TG-AST140041
module load gcc
module load gsl
module list
pwd
date
./mergertree catalogs/catalogNUM0.txt 1 seed0
./mergertree catalogs/catalogNUM1.txt 1 seed1
./mergertree catalogs/catalogNUM2.txt 1 seed2
./mergertree catalogs/catalogNUM3.txt 1 seed3
./mergertree catalogs/catalogNUM4.txt 1 seed4
./mergertree catalogs/catalogNUM5.txt 1 seed5
./mergertree catalogs/catalogNUM6.txt 1 seed6
./mergertree catalogs/catalogNUM7.txt 1 seed7
./mergertree catalogs/catalogNUM8.txt 1 seed8
./mergertree catalogs/catalogNUM9.txt 1 seed9
date
