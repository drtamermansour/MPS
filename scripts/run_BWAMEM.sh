#!/bin/sh

if [ $# -lt 5 ]
then
printf "\nUsage run_BWAMEM.sh [R1 sample] [library Name] [indexed reference] [read path] [script]\n"
exit 0
fi


sample="$1"
lib="$2"
Bwa_ref="$3"
trimmed_reads="$4"
script="$5"

name=$(basename $sample)
SM=$(basename $(dirname $sample))
LB=$SM.$lib
PL="Illumina"                                          
RGID=$(echo $name | cut -d "_" -f1)  ##lane
PU=$RGID.$LB

R1=$trimmed_reads/${name%.fastq.gz}
R2=$(echo $R1 | sed 's/_R1_/_R2_/')
mkdir -p bwa_$RGID
cd bwa_$RGID
echo $Bwa_ref $R1 $R2
echo RGID $RGID LB $LB PL $PL SM $SM
qsub -v RGID="$RGID",SM="$SM",PL="$PL",LB="$LB",PU="$PU",Bwa_ref="$Bwa_ref",R1="$R1",R2="$R2" "${script}"


