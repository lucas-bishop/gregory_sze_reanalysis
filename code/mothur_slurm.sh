#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=lucas_mothur

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=12:00:00

# Account
#SBATCH --account=pschloss1
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=lubi@med.umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=code/log/%x-%j.out

# Environment
##SBATCH --export=ALL

# List compute nodes allocated to the job
if [[ $SLURM_JOB_NODELIST ]] ; then
	echo "Running on"
	scontrol show hostnames $SLURM_JOB_NODELIST
	echo -e "\n"
fi



#####################
#                   #
#  2) Job Commands  #
#                   #
#####################

# Initiating snakemake and running workflow in cluster mode
bash code/mothur/mothurShared.sh data/raw_fastqs/ data/mothur_output/ data/reference_files/
