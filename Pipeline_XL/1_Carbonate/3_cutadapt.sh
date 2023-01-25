#!/bin/bash

#$ -q rcc-30d

# Moger-Reischer_Lennon_XL
#this script uses Slurm. TORQUE/Moab are disabled on Carbonate.
date
cd /N/project/XevoLong
AR=( $(seq 1 88 ) )

for i in "${AR[@]}"
do
    module unload java; module load gatk/3.8; module load samtools; module unload python; module load cutadapt; module load perl/5.30.1; module load fastqc; module load bwa; module unload java; module load picard
	cd Sample_${i}
	
	# clean and trim
	echo "#!/bin/bash" > XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH -J Sample_${i}_qc" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH -p general" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH -o Sample_${i}_trim_stdout.txt" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH -e Sample_${i}_trim_error.txt" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH --mail-type=ALL" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH --mail-user=rzmogerr@indiana.edu" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH --nodes=1" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH --ntasks-per-node=1" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH --cpus-per-task=12" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH --time=2:00:00" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "#SBATCH --mem=20G" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa; module load java; module load picard" >> XL_${i}_clean.sh
	echo ""
	echo "cd /N/project/XevoLong/Sample_${i}" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "time cutadapt -u 17 -u -18 -o Sample_${i}_R1_trimmed.fastq Sample_${i}_R1.fastq" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "time cutadapt -u 17 -u -18 -o Sample_${i}_R2_trimmed.fastq Sample_${i}_R2.fastq" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "mkdir fastqc" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "fastqc -o fastqc/ Sample_${i}_R1_trimmed.fastq" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "fastqc -o fastqc/ Sample_${i}_R2_trimmed.fastq" >> XL_${i}_clean.sh
	echo "" >> XL_${i}_clean.sh
	echo "exit" >> XL_${i}_clean.sh
	
		
	chmod u+x XL_${i}_clean.sh
		
	sbatch XL_${i}_clean.sh
	cd ..
	done