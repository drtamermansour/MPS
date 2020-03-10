#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=24Gb		## 2h x 2 Gb
#mdiag -A ged	
#PBS -m abe			#send email to myself
#PBS -N T_Trim		#give name to the job


module load Trimmomatic/0.33

cd $PBS_O_WORKDIR

#cp $R1_INPUT /tmp/R1.fq
#cp $R2_INPUT /tmp/R2.fq
cd /tmp

java -jar $TRIM/trimmomatic PE -threads 1 -phred33 ${R1_INPUT} ${R2_INPUT} ${output_pe1} ${output_se1} ${output_pe2} ${output_se2} ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:20

cp ${output_pe1} ${output_se1} ${output_pe2} ${output_se2} $PBS_O_WORKDIR/.

qstat -f ${PBS_JOBID}

