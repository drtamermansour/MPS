#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=1:ppn=4,mem=64Gb
#mdiag -A ged
#PBS -m abe
#PBS -N sort_merge

module load SAMTools/1.3.1

cd /tmp

samtools view -u $PBS_O_WORKDIR/pe_aligned_reads.sam | samtools sort -@ 4 - -m 15G -o pe_aligned_reads.sorted.bam
samtools view -u $PBS_O_WORKDIR/se_aligned_reads.sam | samtools sort -@ 4 - -m 15G -o se_aligned_reads.sorted.bam

cd $PBS_O_WORKDIR

rsync -rlpgoDvh --append-verify /tmp/pe_aligned_reads.sorted.bam .
rsync -rlpgoDvh --append-verify /tmp/se_aligned_reads.sorted.bam .


qstat -f ${PBS_JOBID}

