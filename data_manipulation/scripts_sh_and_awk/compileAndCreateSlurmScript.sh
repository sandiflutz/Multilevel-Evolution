#!/bin/bash

if [ -n "$1" ]
then
	macro=$1
	if [ -n "$2" ]
	then
		exename=$2
	else
		exename=exe
	fi	
	if [ -n "$3" ]
	then
		nodo=$3
	else
		nodo=0
	fi	
	script=${exename}.sh
	make cleanobj
	make MACRO=${macro} EXEC=${exename}.out
	echo "#!/bin/bash" > "${script}"
	echo "#SBATCH --job-name=${exename}" >> "${script}"
	echo "#SBATCH --output=output_${exename}.txt" >> "${script}"
	echo "#SBATCH --error=error_${exename}.txt" >> "${script}"
	echo "#SBATCH --time=7-00:00:00" >> "${script}"
	echo "#SBATCH --ntasks=1" >> "${script}"
	echo "#SBATCH --cpus-per-task=1" >> "${script}"
	echo "#SBATCH --partition=${nodo}" >> "${script}"
	echo "cd $SLURM_SUBMIT_DIR" >> "${script}"
	echo "./${exename}.out" >> "${script}"

	echo ""
	echo "To submit the job to the queue, run:"
	echo "		sbatch ${exename}.sh"
fi
