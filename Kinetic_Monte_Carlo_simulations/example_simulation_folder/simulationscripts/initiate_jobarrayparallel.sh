#!/bin/bash
#SBATCH --array=0-600%160
#SBATCH --ntasks=1
#SBATCH --time=80:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --error=error_output_%j.txt
#SBATCH --qos=std
#SBATCH --mail-type=NONE
#SBATCH --mail-user=martijn.vangalen@wur.nl

module load slurm
module load python

inputfile="simulationslist.txt"

declare counter
counter=0
cat $1 | while read basepairs stress repeatnr
do

if [ $counter == $SLURM_ARRAY_TASK_ID ]; then
	
	srun python run_simulation.py ../simulation/BPs_$basepairs/stress_$stress/rep_$repeatnr/ input.txt

fi

counter=$(expr $counter + 1)

done <"$inputfile"
