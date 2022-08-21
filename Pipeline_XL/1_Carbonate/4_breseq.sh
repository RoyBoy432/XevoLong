#!/bin/bash

#$ -q rcc-30d

# Moger-Reischer_Lennon_XL
#this script uses Slurm. TORQUE/Moab are disabled on Carbonate.
date
cd /N/project/XevoLong/breseq
AR=( $(seq 1 79 ) )

for i in "${AR[@]}"
do
    module load breseq
	cd .
	
	# clean and trim
	echo "#!/bin/bash" > XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH -J Sample_${i}_breseq" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH -p general" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH -o Sample_${i}_breseq_stdout.txt" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH -e Sample_${i}_breseq_error.txt" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH --mail-type=ALL" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH --mail-user=rzmogerr@indiana.edu" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH --nodes=1" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH --ntasks-per-node=1" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH --cpus-per-task=12" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH --time=38:00:00" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "#SBATCH --mem=140G" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "module load breseq" >> XL_${i}_breseq.sh
	echo ""
	echo "cd /N/project/XevoLong/breseq" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "breseq -j 8 -p -o breseq_${i} -r /N/project/XevoLong/S288C_NCBI.gb ../Sample_${i}/Sample_${i}_R1_trimmed.fastq ../Sample_${i}/Sample_${i}_R2_trimmed.fastq" >> XL_${i}_breseq.sh
	echo "" >> XL_${i}_breseq.sh
	echo "exit" >> XL_${i}_breseq.sh
	
		
	chmod u+x XL_${i}_breseq.sh
		
	sbatch XL_${i}_breseq.sh
	#cd ..
	done