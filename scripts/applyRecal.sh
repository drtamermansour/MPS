#!/bin/bash -login
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N applyRecal

module load GATK/3.5.0

cd $PBS_O_WORKDIR

java -Xmx45g -jar $GATK/GenomeAnalysisTK.jar \
-T PrintReads \
-R "$gatk_ref" \
-I "$inputBam" \
-BQSR recal_data.table \
-o recal_$inputBam

qstat -f ${PBS_JOBID}

