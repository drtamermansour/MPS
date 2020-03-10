#!/bin/bash -login
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=08:00:00
#SBATCH --mem=48G
#SBATCH -A ged
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH -J GenotypeGVCFs

cd $SLURM_SUBMIT_DIR

module load GATK/4.1.4.1-Python-3.6.4

gatk --java-options "-Xmx40g" GenotypeGVCFs \
-R $gatk_ref \
-V ${sample%.bam}.gatk4.g.vcf \
-O ${sample%.bam}.gatk4.vcf

squeue -l --job ${SLURM_JOB_ID}

