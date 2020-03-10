#!/bin/sh

if [ $# -lt 5 ]
then
printf "\nUsage run_applyRecal.sh [indexed reference fasta] [inputBam] [knownSNPs] [knownINDELs] [script]\n"
exit 0
fi

gatk_ref="$1"
inputBam="$2"
knownSNPs="$3"
knownINDELs="$4"
script="$5"

qsub -v gatk_ref="${gatk_ref}",inputBam="${inputBam}",knownSNPs="${knownSNPs}",knownINDELs="${knownINDELs}" "${script}"
