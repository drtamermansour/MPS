#!/bin/bash -login
#PBS -l walltime=96:00:00,nodes=1:ppn=4,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N haplotypeCaller

module load GATK/3.5.0

cd $PBS_O_WORKDIR

java -Xmx45g -jar $GATK/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $gatk_ref \
$(echo $samples) \
--dbsnp $snps \
-nct 3 \
-o $output
#-o /tmp/$output

#cp /tmp/${output}* $PBS_O_WORKDIR/.

qstat -f ${PBS_JOBID}
