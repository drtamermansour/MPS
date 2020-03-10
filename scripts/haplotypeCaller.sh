#!/bin/bash -login
#PBS -l walltime=48:00:00,nodes=1:ppn=4,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N haplotypeCaller

module load GATK/3.5.0

cd $PBS_O_WORKDIR

java -Xmx45g -jar $GATK/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $gatk_ref \
-I $sample \
-nct 3 \
-o $output ##samplename.vcf ##HC_output_ploidy1_haplo1.vcf

qstat -f ${PBS_JOBID}
