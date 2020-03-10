#!/bin/bash -login
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH -A ged
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH -J haplotypeCaller

cd $SLURM_SUBMIT_DIR

module load GATK/4.1.4.1-Python-3.6.4

gatk --java-options "-Xmx40g" HaplotypeCaller  \
-R $gatk_ref \
-I $sample \
-O ${sample%.bam}.gatk4.g.vcf \
-bamout ${sample%.bam}.bamout.bam \
-ERC GVCF 

squeue -l --job ${SLURM_JOB_ID}

