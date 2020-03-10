#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_genotypeGVCF.sh [known SNPs] [indexed reference fasta] [samples list] [script]\n"
exit 0
fi

snps="$1"
gatk_ref="$2"
sample_list="$3"
output="$4"
script="$5"

samples=""
while read sample; do
  echo $sample
  samples=" ${samples} -I ${sample}"
done < $sample_list

trim_samples=$(echo $samples | xargs | sed 's/\n//')
qsub -v snps="$snps",gatk_ref="${gatk_ref}",samples="${trim_samples}",output="${output}" "${script}"
