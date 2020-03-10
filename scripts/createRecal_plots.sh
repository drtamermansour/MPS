#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N createRecal

module load GATK/3.5.0

cd $PBS_O_WORKDIR

java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R "$gatk_ref" \
-before recal_data.table \
-after post_recal_data.table \
-plots recalibration_plots.pdf

qstat -f ${PBS_JOBID}

