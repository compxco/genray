#!/bin/bash -l 

#SBATCH -J GENRAY
#SBATCH -p regular
#SBATCH --qos=premium
#SBATCH --account=m77
#SBATCH -t 0:02:00
#SBATCH -C haswell
#SBATCH -N 2

## Note: this batchscript can be used at Edison and Cori.
## For Edison, comment the line with haswell.
## Examples of using 480 cores:
##Edison has 24 cores per compute node, 20*24=480, so use (-N 20, -n 160 -c 3)
##Cori has 32 cores per node, 15*32=480, so use (-N 15,  -n 160 -c 3) 

cd $SLURM_SUBMIT_DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/homes/y/ypetrov/pgplot.intel

srun -n 32 -c 4 --cpu_bind=cores ./xgenray_mpi_intel.cori
