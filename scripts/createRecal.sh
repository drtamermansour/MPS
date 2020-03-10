#!/bin/bash -login
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N createRecal

module load GATK/3.5.0

cd $PBS_O_WORKDIR

java -Xmx45g -jar $GATK/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R "$gatk_ref" \
-I "$inputBam" \
-knownSites "$knownSNPs" \
-knownSites "$knownINDELs" \
-o recal_data.table

qstat -f ${PBS_JOBID}

