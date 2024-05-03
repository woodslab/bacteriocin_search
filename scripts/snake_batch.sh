#!/bin/bash

##################
####  Slurm preamble

#### #### ####  These are the most frequently changing options

####  Job name
#SBATCH --job-name=assembly_blast

####  Request resources here
####    These are typically, number of processors, amount of memory,
####    an the amount of time a job requires.  May include processor
####    type, too.

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=13-00:00:00
#SBATCH --account=robertwo1
#SBATCH --output=slurm_logs/%x-%j.log


####  Slurm account and partition specification here
####    These will change if you work on multiple projects, or need
####    special hardware, like large memory nodes or GPUs.

#SBATCH --partition=standard

#### #### ####  These are the least frequently changing options

####  Your e-mail address and when you want e-mail

#SBATCH --mail-user=agarret@umich.edu
#SBATCH --mail-type=BEGIN,FAIL,END

####  End Slurm preamble


##################
#conda activate new
#xport OMP_NUM_THREADS=24
#echo $OMP_NUM_THREADS
#conda activate snakemake

snakemake -s Snakefile -j 96 --use-conda --rerun-incomplete



 
