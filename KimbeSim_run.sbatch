#!/bin/bash -l
#SBATCH --time=23:00:00
#SBATCH --nodes=1 # stick with one node; special instructions (MPI) are required to communicate among nodes
#SBATCH --cpus-per-task=1 # i.e., cores per task
#SBATCH --ntasks=1
#SBATCH --mem=275g
#SBATCH --output experiments/Exp4_20260611/output/slurm_out/%A_%a.output
#SBATCH --error experiments/Exp4_20260611/output/slurm_out/%A_%a.output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schla103@umn.edu
# #SBATCH --array=1   # task id is ${SLURM_ARRAY_TASK_ID}
# #SBATCH --dependency=afterok:11092220

module load geos/3.11.2-rocky8
module load gdal/3.8.2-rocky8
module load proj/9.3.1-rocky8
module load udunits/2.2.27.6_gcc11.3.0-rocky8
module load openssl/3.1.1
module load R/4.4.0-openblas-rocky8

Rscript KimbeSim_run.R
