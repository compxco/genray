#!/bin/bash -l 

#Invoke as: sbatch edison_batchscript_mpi.slurm
#Check run status with: sqs


#SBATCH -J GENRAY
#SBATCH -p debug
#SBATCH --account=m77
#SBATCH -t 0:05:00
#SBATCH -N 1
#Edison has 24 cores per compute node , so 1*24=24

#test7 launches 5 rays.  Parallelize with number of
#  rays plus 1 (controller) core. 

srun -n 6  /global/homes/u/u650/genray_cswim_svn/trunk/xgenray_mpi_intel.edison > log_test7_5rays

#BH180907
#This run takes full time: 22.45 secs.   To be compared with serial version
#  which takes 58.2 secs.   Since parallelization is over rays, then the
#  execution time for the mpi version is controlled by the ray which takes
#  the longest time.
