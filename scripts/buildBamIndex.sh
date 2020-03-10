#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=2,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N buildBamIndex

module load picardTools/1.113

cd $PBS_O_WORKDIR

java -Xmx45g -jar $PICARD/BuildBamIndex.jar INPUT=$sample

qstat -f ${PBS_JOBID}
