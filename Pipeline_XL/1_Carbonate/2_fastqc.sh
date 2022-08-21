#!/bin/bash

#$ -q rcc-30d

# Moger-Reischer_Lennon_XL
#this script uses Slurm. TORQUE/Moab are disabled on Carbonate.
date
cd /N/project/XevoLong
AR=( $(seq 1 88 ) )

for i in "${AR[@]}"
do
    module load gatk/3.8; module load samtools; module unload python; module load cutadapt; module load perl/5.30.1; module load fastqc; module load bwa; module unload java; module load picard
	cd Sample_${i}
	cat *R1_001.fastq > Sample_${i}_R1.fastq
	cat *R2_001.fastq > Sample_${i}_R2.fastq
	
	# run qc
	echo "#!/bin/bash" > XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH -J Sample_${i}_qc" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH -p general" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH -o Sample_${i}_qc_stdout.txt" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH -e Sample_${i}_qc_error.txt" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH --mail-type=ALL" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH --mail-user=rzmogerr@indiana.edu" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH --nodes=1" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH --ntasks-per-node=1" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH --cpus-per-task=12" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH --time=2:00:00" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "#SBATCH --mem=20G" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "module load samtools; module unload python; module load python/2.7.16; module load cutadapt; module load perl/5.30.1; module load fastqc; module load bwa; module load java; module load picard" >> XL_${i}_qc.sh
	echo ""
	echo "cd /N/project/XevoLong/Sample_${i}" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "mkdir fastqc_untrimmed" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "fastqc -o fastqc_untrimmed/ Sample_${i}_R1.fastq" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "fastqc -o fastqc_untrimmed/ Sample_${i}_R2.fastq" >> XL_${i}_qc.sh
	echo "" >> XL_${i}_qc.sh
	echo "exit" >> XL_${i}_qc.sh
	
		
	chmod u+x XL_${i}_qc.sh

		
	sbatch XL_${i}_qc.sh
	cd ..
	done