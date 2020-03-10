#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_snpVariantFiltration.sh [indexed reference fasta] [inputVcf] [script]\n"
exit 0
fi

gatk_ref="$1"
inputVcf="$2"
outputVcf="$3"
script="$4"

qsub -v gatk_ref="${gatk_ref}",inputVcf="${inputVcf}",outputVcf="${outputVcf}" "${script}"




