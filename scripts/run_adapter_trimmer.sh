#!/bin/sh

if [ $# -lt 2 ]
then
    printf "\nUsage run_adapter_trimmer.sh [R1 Sample] [script]\n"
    exit 0
fi

f="$1"
script="$2"

echo $f
input_one="$f"
input_two=$(echo "$f" | sed s/_R1_/_R2_/)
temp=$(basename "$f")
temp2=${temp%.fastq.gz}
output_pe1=$temp2".pe.fq"
output_pe2=$(echo "$output_pe1" | sed s/_R1_/_R2_/)
output_se1=$temp2".se.fq"
output_se2=$(echo "$output_se1" | sed s/_R1_/_R2_/)
qsub -v R1_INPUT="$input_one",R2_INPUT="$input_two",output_pe1="$output_pe1",output_pe2="$output_pe2",output_se1="$output_se1",output_se2="$output_se2" $script





