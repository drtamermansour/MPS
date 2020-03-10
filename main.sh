## define pathes
work_dir="/mnt/home/mansourt/Tamer/MPS"
mkdir -p $work_dir/{scripts,fastq_data,varResults,varResults_snps,varResults_indels}
script_path=$work_dir/scripts
varResults=$work_dir/varResults
varResults_snps=$work_dir/varResults_snps
varResults_indels=$work_dir/varResults_indels

dogSeq="/mnt/home/mansourt/Tamer/dogSeq"
genome_dir=$dogSeq/refGenome

dogAnn="/mnt/home/mansourt/Tamer/dogAnn"
dogAnnRes=$dogAnn/refResources

mkdir -p $work_dir/fastq_data/Ellie
cd $work_dir/fastq_data/Ellie

wget http://slims.bioinformatics.ucdavis.edu/Data/hbs49qnt6/Unaligned/Project_TMTM_L8_Ellie_MPS/Ellie-MPS_S127_L008_R1_001.fastq.gz
wget http://slims.bioinformatics.ucdavis.edu/Data/hbs49qnt6/Unaligned/Project_TMTM_L8_Ellie_MPS/Ellie-MPS_S127_L008_R2_001.fastq.gz
wget http://slims.bioinformatics.ucdavis.edu/Data/hbs49qnt6/Unaligned/Project_TMTM_L8_Ellie_MPS/Ellie-MPS_S127_L008_R3_001.fastq.gz
wget http://slims.bioinformatics.ucdavis.edu/Data/hbs49qnt6/Unaligned/Project_TMTM_L8_Ellie_MPS/laneBarcode.html
vim runInfo.txt ## the view run page data
mkdir indexData
mv Ellie-MPS_S127_L008_R2_001.fastq.gz indexData/Ellie-MPS_S127_L008_index.fastq.gz
mv Ellie-MPS_S127_L008_R3_001.fastq.gz Ellie-MPS_S127_L008_R2_001.fastq.gz
chmod u-w *

## fastqc
bash ${script_path}/run_fastqc.sh "$work_dir/fastq_data";

## Adapter trimming
## Mild trimming with Trimmomatic using sliding window
R1=$work_dir/fastq_data/Ellie/Ellie-MPS_S127_L008_R1_001.fastq.gz
mkdir -p $work_dir/trimmed_RNA_reads
cd $work_dir/trimmed_RNA_reads
bash ${script_path}/run_adapter_trimmer.sh $R1 $script_path/adapter_trimmer_PE.sh

## get the referenece genome
#cd $genome_dir
#wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz' -O canFam3.fa.gz
#gunzip canFam3.fa.gz

## prepare BWA index (for Mapping reads)
#mkdir -p $genome_dir/BwaIndex && cd $genome_dir/BwaIndex
#cp ../canFam3.fa genome.fa
#bash ${script_path}/run_bwa-index.sh genome.fa
Bwa_ref="$genome_dir/BwaIndex/genome.fa"

## read mapping
mkdir -p $work_dir/bwa_align
cd $work_dir/bwa_align
bash ${script_path}/run_BWAMEM.sh "$R1" "lib1" "$Bwa_ref" "$work_dir/trimmed_RNA_reads" "${script_path}/BWAMEM.sh"
## Convert to BAM and sort
cd bwa_Ellie-MPS
qsub ${script_path}/getBAMfile.sh
## build bam index
sample="pe_aligned_reads.sorted.bam"
qsub -v sample=${sample} $script_path/buildBamIndex.sh
## assess sample coverage 
module load SAMTools/1.3.1
samtools depth $sample | awk '{sum+=$3} END { print "Average = ",sum/NR}' > sampleCoverage
## assess mapping rate
samtools flagstat pe_aligned_reads.sorted.bam
## assess insert size
module load picardTools/1.89
java -jar $PICARD/CollectInsertSizeMetrics.jar I=pe_aligned_reads.sorted.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf 

## prepare GATK dictionary and index (for GATK variant analysis)
#mkdir -p $genome_dir/gatkIndex && cd $genome_dir/gatkIndex
#cp ../canFam3.fa genome.fa
#bash ${script_path}/run_gatk-index.sh genome.fa
gatk_ref="$genome_dir/gatkIndex/genome.fa"
gatk_ref_index="$genome_dir/gatkIndex/genome.fa.fai"

## create known var
dogSeqSNPs=$dogSeq/varResults/GenotypeGVCFs_output_max50.raw_SNPs.vcf
dogSeqINDELs=$dogSeq/varResults/GenotypeGVCFs_output_max50.raw_INDELs.vcf
mkdir $work_dir/knownVar && cd $work_dir/knownVar
bash ${script_path}/run_snpVariantFiltration.sh "$gatk_ref" "$dogSeqSNPs" "hundredDogGenome.filtred_snps.vcf" "$script_path/snpVariantFiltration_aggressive.sh"
bash ${script_path}/run_indelVariantFiltration.sh "$gatk_ref" "$dogSeqINDELs" "hundredDogGenome.filtred_indels.vcf" "$script_path/indelVariantFiltration_aggressive.sh"

#module load vcftools/0.1.14
#vcftools --vcf hundredDogGenome.filtred_snps.vcf --remove-filtered-all --recode --out hundredDogGenome_snps --max-indv 0
module load bcftools/1.2
bcftools view --drop-genotypes -f .,PASS -Ov -o hundredDogGenome_snps2.vcf hundredDogGenome.filtred_snps.vcf ## 11034217
bcftools view --drop-genotypes -f .,PASS -Ov -o hundredDogGenome_indels2.vcf hundredDogGenome.filtred_indels.vcf ## 4947081
knownSNPs=$work_dir/knownVar/hundredDogGenome_snps2.vcf
knownINDELs=$work_dir/knownVar/hundredDogGenome_indels2.vcf

## base recalibration
inputBam="pe_aligned_reads.sorted.bam"
cd $work_dir/bwa_align/bwa_Ellie-MPS
bash ${script_path}/run_createRecal.sh "$gatk_ref" "$inputBam" "$knownSNPs" "$knownINDELs" "$script_path/createRecal.sh"
bash ${script_path}/run_createRecal.sh "$gatk_ref" "$inputBam" "$knownSNPs" "$knownINDELs" "$script_path/createRecal_2nd.sh"
bash ${script_path}/run_createRecal.sh "$gatk_ref" "$inputBam" "$knownSNPs" "$knownINDELs" "$script_path/createRecal_plots.sh"
bash ${script_path}/run_createRecal.sh "$gatk_ref" "$inputBam" "$knownSNPs" "$knownINDELs" "$script_path/applyRecal.sh"

## variant calling
sample="recal_pe_aligned_reads.sorted.bam"
output=${sample%.bam}.vcf
qsub -v gatk_ref="${gatk_ref}",sample="${sample}",output="${output}" $script_path/haplotypeCaller.sh
module load bcftools/1.2
bcftools convert -Oz -o recal_pe_aligned_reads.sorted.vcf.gz recal_pe_aligned_reads.sorted.vcf
bcftools index --tbi recal_pe_aligned_reads.sorted.vcf.gz

## joint variant calling with control genomes
## define the list samples.
#echo $work_dir/bwa_align/bwa_Ellie-MPS/recal_pe_aligned_reads.sorted.bam > MPSwithControl.list
#echo $dogSeq/data/newSeq/bwa_align/bwa_3726/dedup_reads.bam >> MPSwithControl.list
#echo $dogSeq/data/newSeq/bwa_align/bwa_4885/dedup_reads.bam >> MPSwithControl.list
## joint genotyping
#sample_list="MPSwithControl.list"
#output="MPSwithControl_jointCalling.vcf"
#bash ${script_path}/run_haplotypeCaller.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$output" "$script_path/haplotypeCaller_multi.sh"

## get control VCF
dogSeqVar=$dogSeq/varResults/GenotypeGVCFs_output_max50.vcf
cd $work_dir/knownVar
module load bcftools/1.2
## control VCF for using the 2 Boston_Terrier dogs
bcftools view --samples 3726,4885 --min-ac 1 -Ov -o control_Boston_Terrier.vcf $dogSeqVar
bcftools convert -Oz -o control_Boston_Terrier.vcf.gz control_Boston_Terrier.vcf
bcftools index --tbi control_Boston_Terrier.vcf.gz
## control VCF for all the 100 genome dogs 
bcftools convert -Oz -o control_dogSeq.vcf.gz $dogSeqVar
bcftools index --tbi control_dogSeq.vcf.gz

## merge with Boston_Terrier control files (we ended up not using these files)
# first approach: using bcftools
cd $work_dir/bwa_align/bwa_Ellie-MPS
bcftools merge --merge all -Ov -o MPSwithControl.vcf $work_dir/knownVar/control_Boston_Terrier.vcf.gz recal_pe_aligned_reads.sorted.vcf.gz

# second approach: using GATK
cd $work_dir/bwa_align/bwa_Ellie-MPS
module load GATK/3.5.0
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $gatk_ref \
--variant:ensembl $work_dir/knownVar/control_Boston_Terrier.vcf \
--variant:broad recal_pe_aligned_reads.sorted.vcf \
-o MPSwithControl2.vcf \
-genotypeMergeOptions UNIQUIFY

## compare files
module load vcftools/0.1.14
#module load zlib/1.2.8
vcftools --vcf $dogSeqVar --gzdiff recal_pe_aligned_reads.sorted.vcf.gz --diff-site --out allControls_v_MPS
#vcftools --gzvcf $work_dir/knownVar/control_Boston_Terrier.vcf.gz --gzdiff input_file2.vcf.gz --diff-site --out in1_v_in2

## disease specific varaint annotation
broad_gtf=$dogAnnRes/canis_familiaris.protein_coding.gtf
ens_gtf=$dogAnnRes/Canis_familiaris.CanFam3.1.86.gtf
ncbi_gtf=$dogAnnRes/ref_CanFam3.1_top_level_mapped.gtf
ncbi_gff=$dogAnnRes/ref_CanFam3.1_top_level.gff3
ncbiMapped_gff=$dogAnnRes/ref_CanFam3.1_top_level_mapped.gff3
genome_ens=$genome_dir/canFam3_ens.fa
###########
## candidate gene approach 
## MPS
##get gene co-ordinates
label="MPS"  ## label="ML"
mkdir -p $work_dir/${label}_genes && cd $work_dir/${label}_genes ## create file MPS_genes.txt from litrature 

#while read gene;do grep $gene $ncbi_gff >> ${label}.refGene.gff3;done < ${label}_genes.txt
#cat ${label}.refGene.gff3 | awk -F "\t" '$3=="mRNA"' | grep "Note=" > ${label}.Note.gff3 ## 3 transcripts changed compared to genomic sequence

#mkdir geneSp
#while read gene;do grep $gene $broad_gtf > geneSp/$gene.broad.gtf;done < ${label}_genes.txt
#while read gene;do grep $gene $ncbi_gtf > geneSp/$gene.refGene.gtf;done < ${label}_genes.txt
#while read gene;do grep $gene $ncbi_gtf > geneSp/$gene.ensembl.gtf;done < ${label}_genes.txt
#> ${label}.broad.gtf; while read gene;do grep $gene $broad_gtf >> ${label}.broad.gtf;done < ${label}_genes.txt;
#egrep -v "GLB1L|FARSB|DAGLB|SH3GLB" ${label}.broad.gtf > ${label}.broad_v2.gtf
> ${label}.refGene.gtf; while read gene;do grep $gene $ncbi_gtf >> ${label}.refGene.gtf;done < ${label}_genes.txt;
egrep -v "GLB1L|FARSB|DAGLB|SH3GLB" ${label}.refGene.gtf > ${label}.refGene_v2.gtf
> ${label}.ensembl.gtf; while read gene;do grep $gene $ens_gtf >> ${label}.ensembl.gtf;done < ${label}_genes.txt;
egrep -v "GLB1L|FARSB|DAGLB|SH3GLB" ${label}.ensembl.gtf > ${label}.ensembl_v2.gtf
for f in ${label}.*_v2.gtf;do bash $script_path/gtfToBed.sh "$f" "$script_path";done
awk -F "\t" -v OFS='\t' '{ print "chr"$0 }' ${label}.ensembl_v2.bed > ${label}.ensembl2_v2.bed
#cat ${label}.broad_v2.bed ${label}.refGene_v2.bed ${label}.ensembl2_v2.bed  > all_${label}.bed
cat ${label}.refGene_v2.bed ${label}.ensembl2_v2.bed  > all_${label}.bed
sort -k1,1 -k2,2n all_${label}.bed > all_${label}.sorted.bed
module load BEDTools/2.24.0
bedtools merge -i all_${label}.sorted.bed > all_${label}.reduced.bed

## prepare local caches for off-line variant annotation
module load VEP/85
gtf2vep.pl -i ${label}.ensembl_v2.gtf -f $genome_ens -d 85 -s canFam_${label}.Ensembl --verbose &> gtf2vep_ens.log

#cat ${label}.refGene_v2.gtf | awk -F "\t" -v OFS='\t' '{ print $0" gene_source \"NCBI\"; gene_biotype \"protein_coding\"; transcript_source \"NCBI\"; transcript_biotype \"protein_coding\";" }' > ${label}.refGene_VEP.gtf
#while read line;do
# echo "$line" | awk -F "\t" -v OFS='\t' '{if($3=="start_codon" || $3=="stop_codon")print;}' | sed 's/exon_id ".*"; gene_name/gene_name/';
# echo "$line" | awk -F "\t" -v OFS='\t' '{if($3!="start_codon" && $3!="stop_codon")print;}'
#done < ${label}.refGene_VEP.gtf > ${label}.refGene_VEP2.gtf
#sed 's/^chr//' ${label}.refGene_VEP2.gtf > ${label}.refGene_VEP3.gtf
#gtf2vep.pl -i ${label}.refGene_VEP3.gtf -f $genome_ens -d 85 -s canFam_${label}.NCBI --no_transcripts --verbose &> gtf2vep_ncbi.log
> ${label}.refGene.gff; while read gene;do grep $gene $ncbiMapped_gff >> ${label}.refGene.gff;done < ${label}_genes.txt;
egrep -v "GLB1L|FARSB|DAGLB|SH3GLB" ${label}.refGene.gff > ${label}.refGene_v2.gff
sed 's/^chr//' ${label}.refGene_v2.gff > ${label}.refGene_VEP.gff
gtf2vep.pl -i ${label}.refGene_VEP.gff -f $genome_ens -d 85 -s canFam_${label}.NCBIgff --verbose &> gtf2vep_ncbiGFF.log         
#gtf2vep.pl -i ${label}.refGene_VEP.gff -f $genome_ens -d 85 -s canFam_${label}.NCBIgff --no_transcripts --verbose &> gtf2vep_ncbiGFF.log

#cat ${label}.broad_v2.gtf | awk -F "\t" -v OFS='\t' '{ print $0" gene_source \"Broad\"; gene_biotype \"protein_coding\"; transcript_source \"Broad\"; transcript_biotype \"protein_coding\";" }' > ${label}.broad_VEP.gtf
#while read line;do
# echo "$line" | awk -F "\t" -v OFS='\t' '{if($3=="start_codon" || $3=="stop_codon")print;}' | sed 's/exon_id ".*"; gene_source/gene_source/';
# echo "$line" | awk -F "\t" -v OFS='\t' '{if($3!="start_codon" && $3!="stop_codon")print;}'
#done < ${label}.broad_VEP.gtf > ${label}.broad_VEP2.gtf
#sed 's/^chr//' ${label}.broad_VEP2.gtf > ${label}.broad_VEP3.gtf
#gtf2vep.pl -i ${label}.broad_VEP3.gtf -f $genome_ens -d 85 -s canFam_${label}.Broad --no_transcripts --verbose &> gtf2vep_broad.log

cd $work_dir/bwa_align/bwa_Ellie-MPS
vcfFile_gz="recal_pe_aligned_reads.sorted.vcf.gz"
module load bcftools/1.2
bcftools view -Ov -R $work_dir/${label}_genes/all_${label}.reduced.bed -o ${label}_specific.vcf $vcfFile_gz  ## 968 variants
bcftools view -Ov -R $work_dir/${label}_genes/all_${label}.reduced.bed --min-ac 2 -o ${label}_specific_homo.vcf $vcfFile_gz ## 482 homo varaints
#grep '^#' ${label}_specific_homo.vcf > ${label}_specific_homo_sorted.vcf && grep -v '^#' ${label}_specific_homo.vcf | LC_ALL=C sort -t '\t' -k1,1 -k2,2n >> ${label}_specific_homo_sorted.vcf
variant_effect_predictor.pl -offline -i ${label}_specific_homo.vcf --species canFam_${label}.Ensembl -o ${label}_homo.spEnsembl.ann
variant_effect_predictor.pl -offline -i ${label}_specific_homo.vcf --species canFam_${label}.NCBIgff -o ${label}_homo.spNCBI.ann
#variant_effect_predictor.pl -offline -i ${label}_specific_homo.vcf --species canFam_${label}.Broad -o ${label}_homo.spBroad.ann
for f in *.ann;do grep "IMPACT=HIGH" $f;done > highImpact
awk -v ORF="\t" '{print $1}' highImpact | sort | uniq | wc -l ## only 1

## compare to control_Boston_Terrier.vcf.gz
module load vcftools/0.1.14
bcftools view -Ov -R $work_dir/${label}_genes/all_${label}.reduced.bed -o Boston_specific.vcf $work_dir/knownVar/control_Boston_Terrier.vcf.gz
vcftools --vcf ${label}_specific_homo.vcf --diff Boston_specific.vcf --diff-site --out ${label}_homo.vs.Boston
head -n1 ${label}_homo.vs.Boston.diff.sites_in_files > ${label}_specific_homo.vs.Boston
awk -F "\t" -v OFS='\t' '{if($4==1)print;}' ${label}_homo.vs.Boston.diff.sites_in_files >> ${label}_specific_homo.vs.Boston ## 74 varaint                                                

## comare to control_dogSeq.vcf.gz
bcftools view -Ov -R $work_dir/${label}_genes/all_${label}.reduced.bed -o dogSeq_specific.vcf $work_dir/knownVar/control_dogSeq.vcf.gz
vcftools --vcf ${label}_specific_homo.vcf --diff dogSeq_specific.vcf --diff-site --out ${label}_homo.vs.Allcontrol 
head -n1 ${label}_homo.vs.Allcontrol.diff.sites_in_files > ${label}_specific_homo.vs.all
awk -F "\t" -v OFS='\t' '{if($4==1)print;}' ${label}_homo.vs.Allcontrol.diff.sites_in_files >> ${label}_specific_homo.vs.all ## 7 varaint
###########
#### whole genome approach
cd $work_dir/bwa_align/bwa_Ellie-MPS
## get homozygous variants in Ellie
#bcftools view -Ov --min-ac 2 -o ${label}_homo.vcf $vcfFile_gz ## 2560530 homo varaints

## get homozygous variants in Ellie but exclude uncharacterized chromosomes 
head -n39 $dogAnnRes/ncbi/chromosome_map.txt | awk -v OFS="\t" '{print $5,"1",$6}' > pos.txt
bcftools view -Ov -R pos.txt --min-ac 2 -o ${label}_homo.vcf $vcfFile_gz ## 2498167 homo varaints

## filter by variants in Boston_Terrier
bcftools view -Ov -R pos.txt -o Boston.vcf $work_dir/knownVar/control_Boston_Terrier.vcf.gz   
vcftools --vcf ${label}_homo.vcf --diff Boston.vcf --diff-site --out ${label}_all_homo.vs.Boston
head -n1 ${label}_all_homo.vs.Boston.diff.sites_in_files > ${label}_homo.vs.Boston
awk -F "\t" -v OFS='\t' '{if($4==1)print;}' ${label}_all_homo.vs.Boston.diff.sites_in_files >> ${label}_homo.vs.Boston # 464589 varaint

## filter by variants in 101 genome
bcftools view -Ov -R pos.txt -o dogSeq.vcf $work_dir/knownVar/control_dogSeq.vcf.gz
vcftools --vcf ${label}_homo.vcf --diff dogSeq.vcf --diff-site --out ${label}_all_homo.vs.Allcontrol
head -n1 ${label}_all_homo.vs.Allcontrol.diff.sites_in_files > ${label}_homo.vs.all   
awk -F "\t" -v OFS='\t' '{if($4==1)print;}' ${label}_all_homo.vs.Allcontrol.diff.sites_in_files >> ${label}_homo.vs.all # 31648 varaint
tail -n+2 ${label}_homo.vs.all | awk -F "\t" -v OFS='\t' '{print $1,$2;}' > candidate.list
vcftools --vcf ${label}_homo.vcf --positions candidate.list --recode --out ${label}_homo.vs.all

#### functional annotation of the filtered list
## prepare the VCF_info tables (VCF dependent)
echo -e "Location\tCHROM\tPOS\tID\tREF\tALT" > ${label}_VCF_info
grep -v "^chrUn_" ${label}_homo.vs.all.recode.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{print substr($1,4)":"$2,$1,$2,$3,$4,$5}' >> ${label}_VCF_info

## prepare extended annotation files
extEnsAnn=$dogSeq/varResults/breedSp/ENS_TransInfo.txt
extNCBIAnn=$dogSeq/varResults/breedSp/NCBI_TransInfo.txt

## prepare variant effect files
module load VEP/85
## Ensembl
#variant_effect_predictor.pl -offline -i ${label}_homo.vs.all.recode.vcf --species canFam.Ensemblgff -o ${label}_homo.vs.all.Ensemblgff.ann 
#variant_effect_predictor.pl -offline -i ${label}_homo.vs.all.recode.vcf --species canFam.Ensembl -o ${label}_homo.vs.all.Ensembl.ann
variant_effect_predictor.pl -i ${label}_homo.vs.all.recode.vcf --cache --offline --cache_version 86 -species canis_familiaris -o ${label}_homo.vs.all.Ensembldb.ann

# add the location field
grep -v "^##" ${label}_homo.vs.all.Ensembldb.ann | sed 's/#Uploaded_variation/Uploaded_variation/' > tempVarEff.txt
echo "Location" > tempVarEff2.txt;
grep -v "^#" ${label}_homo.vs.all.Ensembldb.ann | awk 'BEGIN{FS="[\t:]";}{split($3,a,"-");if($4=="-"){print $2":"a[1]-1}else{print $2":"a[1]}}' >> tempVarEff2.txt;
paste tempVarEff2.txt tempVarEff.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' > ${label}_homo.vs.all.Ensembldb.ann2
rm tempVarEff*.txt

# select high impact variants
grep "IMPACT=HIGH" ${label}_homo.vs.all.Ensembldb.ann2 > highImpact.Ensembldb
awk -v ORF="\t" '{print $1}' highImpact.Ensembldb | sort | uniq | wc -l ## 75

# extend the annotation
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE,quote="");'\
'dataMerge=merge(data1,data2,by.x="V5",by.y="Transcript_ID",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"ext",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact.Ensembldb" $extEnsAnn 

# add VCF info
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE);'\
'dataMerge=merge(data1,data2,by.x="V3",by.y="Location",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"VCFinfo",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact.Ensembldb.ext" ${label}_VCF_info

## NCBI
variant_effect_predictor.pl -offline -i ${label}_homo.vs.all.recode.vcf --species canFam.NCBIgff -o ${label}_homo.vs.all.NCBI.ann

# add the location field
grep -v "^##" ${label}_homo.vs.all.NCBI.ann | sed 's/#Uploaded_variation/Uploaded_variation/' > tempVarEff.txt
echo "Location" > tempVarEff2.txt;
grep -v "^#" ${label}_homo.vs.all.NCBI.ann | awk 'BEGIN{FS="[\t:]";}{split($3,a,"-");if($4=="-"){print $2":"a[1]-1}else{print $2":"a[1]}}' >> tempVarEff2.txt;
paste tempVarEff2.txt tempVarEff.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' > ${label}_homo.vs.all.NCBI.ann2
rm tempVarEff*.txt

# select high impact variants
grep "IMPACT=HIGH" ${label}_homo.vs.all.NCBI.ann2 > highImpact.NCBI
awk -v ORF="\t" '{print $1}' highImpact.NCBI | sort | uniq | wc -l ## 49 

# find gene IDs
awk -F "\t" '{print $1,$5,$7}' highImpact.NCBI | sort | uniq > highImpact.NCBI_transcripts  ## the 49 muation affect 39 transcripts
while read var trans eff;do grep $trans $ncbi_gff | head -n1 | awk -F "[\t=;,]" '{print $14}';done < highImpact.NCBI_transcripts > highImpact.NCBI_geneIDS ## in 34 gene
paste highImpact.NCBI_transcripts highImpact.NCBI_geneIDS > highImpact.NCBI_genes

# find genes with assembly notes
#while read id;do grep $id $ncbi_gff | grep "Note=" | awk -F "[\t=;,]" '{print $14}' | uniq;done < highImpact.NCBI_geneIDS | uniq > highImpact.NCBI_geneIDS.genomeEdit
#cat highImpact.NCBI_geneIDS.genomeEdit | grep -Fwf - highImpact.NCBI_genes > highImpact.NCBI_genes.genomeEdit  ## 21 variants affecting 12 transcript in 12 genes
#cat highImpact.NCBI_geneIDS.genomeEdit | grep -v -Fwf - highImpact.NCBI_genes > highImpact.NCBI_genes.candidates  ## 28 variants affecting 27 transcript in 22  genes (5 genes has more than 1 mutation suggesting wrong genome sequence or wrong model
while read id;do grep $id $ncbi_gff | grep "Note=" | awk -F "\t" '{print $9}' | sed 's/.*GeneID/GeneID/' | sed 's/,.*;Note/\tNote/' | sed 's/;.*//' | grep "RefSeq protein" | sort | uniq;done < highImpact.NCBI_geneIDS | sort | uniq > highImpact.NCBI_geneIDS.genomeEdit_details
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1]);data2=read.table(args[2],sep="\t");'\
'dataMerge=merge(data1,data2,by.x="V4",by.y="V1",all.x = TRUE);'\
'write.table(unique(dataMerge[,c(1,3,5)]),paste(args[1],"details",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact.NCBI_genes" "highImpact.NCBI_geneIDS.genomeEdit_details" ## merge by gene ID

# merge annotation file with gene assembly notes
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t");'\
'dataMerge=merge(data1,data2,by.x="V5",by.y="V2",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"Notes",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact.NCBI" "highImpact.NCBI_genes.details" ## merge by transcript ID

# extend the annotation
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE,quote="");'\
'dataMerge=merge(data1,data2,by.x="V1",by.y="Transcript_ID",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"ext",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact.NCBI.Notes" $extNCBIAnn

# add VCF info
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE);'\
'dataMerge=merge(data1,data2,by.x="V3",by.y="Location",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"VCFinfo",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact.NCBI.Notes.ext" ${label}_VCF_info

## shared annotation between Ensembl and NCBI
awk -v ORF="\t" '{print $1}' highImpact.Ensembldb | uniq > highImpact_ids.Ensembldb
awk -v ORF="\t" '{print $1}' highImpact.NCBI | uniq > highImpact_ids.NCBI
comm <(sort highImpact_ids.NCBI) <(sort highImpact_ids.Ensembldb) > highImpact.comp  ## 14 are shared


###################################################
######## Re-analysis after the paper review 
module load bcftools/1.9.64  ## module load bcftools/1.2
module load VCFtools/0.1.15-Perl-5.26.1 ## module load vcftools/0.1.14
## get control VCF with a high allele frequency (> 1%) 
dogSeqVar=$dogSeq/varResults/GenotypeGVCFs_output_max50.vcf
cd $work_dir/knownVar
## control VCF for all the 100 genome dogs 
bcftools view --min-af 0.01 -Ov -o control_dogSeq_0.01.vcf $dogSeqVar ## remove 5,429,659 from 24,193,426 variants
bcftools convert -Oz -o control_dogSeq_0.01.vcf.gz control_dogSeq_0.01.vcf
bcftools index --tbi control_dogSeq_0.01.vcf.gz

#### whole genome approach
cd $work_dir/bwa_align/bwa_Ellie-MPS

## get homozygous variants in Ellie but exclude uncharacterized chromosomes (No change)
head -n39 $dogAnnRes/ncbi/chromosome_map.txt | awk -v OFS="\t" '{print $5,"1",$6}' > pos.txt
bcftools view -Ov -R pos.txt --min-ac 2 -o ${label}_homo.vcf $vcfFile_gz ## 2498167 homo varaints
#grep -v "^#" ${label}_homo.vcf | grep -v "AC=2;" 

## filter by variants in 101 genome with high MAF (>1%)
bcftools view -Ov -R pos.txt -o dogSeq_0.01.vcf $work_dir/knownVar/control_dogSeq_0.01.vcf.gz
vcftools --vcf ${label}_homo.vcf --diff dogSeq_0.01.vcf --diff-site --out ${label}_all_homo.vs.Allcontrol_0.01
head -n1 ${label}_all_homo.vs.Allcontrol_0.01.diff.sites_in_files > ${label}_homo.vs.all_0.01
awk -F "\t" -v OFS='\t' '{if($4==1)print;}' ${label}_all_homo.vs.Allcontrol_0.01.diff.sites_in_files >> ${label}_homo.vs.all_0.01 # 34921 instead of 31648 varaint
tail -n+2 ${label}_homo.vs.all_0.01 | awk -F "\t" -v OFS='\t' '{print $1,$2;}' > candidate.list_0.01
vcftools --vcf ${label}_homo.vcf --positions candidate.list_0.01 --recode --out ${label}_homo.vs.all_0.01
#grep -v "^#" ${label}_homo.vs.all_0.01.recode.vcf | grep -v "1/1:"


#### functional annotation of the filtered list
## prepare the VCF_info tables (VCF dependent)
echo -e "Location\tCHROM\tPOS\tID\tREF\tALT" > ${label}_VCF_info_0.01
grep -v "^chrUn_" ${label}_homo.vs.all_0.01.recode.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{print substr($1,4)":"$2,$1,$2,$3,$4,$5}' >> ${label}_VCF_info_0.01

## A) Extend ENSEMBL annotation:
# Download tables for dog genes from Biomart (www.ensembl.org/biomart) >> Database:Ensembl Genes 98, Dataset:Canis familiaris genes(CanFam3.1) ,Attributes > GENE > Ensembl > Ensembl Transcript ID, Associated Gene Name, Gene type, Description. Select results export to "file", with "TSV" format, and select "Unique results only". Save as "dogEnsembl_TransInfo.txt"
echo -e "Transcript_ID\tGene_Name\tGene_biotype\tDescription" > ENS_TransInfo.txt
tail -n+2 dogEnsembl_TransInfo.txt | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "-" }; { print $1,$3,$4,$2 }' >> ENS_TransInfo.txt ## dogEnsembl_TransInfo_noEmpty.txt
## B) Extend NCBI annotation
wget ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/GFF/ref_CanFam3.1_top_level.gff3.gz
gunzip ref_CanFam3.1_top_level.gff3.gz
grep "ID=rna" ref_CanFam3.1_top_level.gff3 | awk -F "\t" '{print $9}' | awk -F "[,;]" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; Name = vars["Name"]; gene = vars["gene"]; product = vars["product"]; } { print Dbxref,Name,gene,product }' > NCBI_TransInfo.temp.trans;
grep -v "^#" ref_CanFam3.1_top_level.gff3 | awk -F "\t" '{if($3=="gene")print $9}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; gene_biotype = vars["gene_biotype"]; } { print Dbxref,gene_biotype }' > NCBI_TransInfo.temp.gene;
Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1],header=F);data2=read.delim(args[2],header=F);'\
'dataMerge=merge(data1,data2,by="V1",all.x=F,all.y=T); colnames(dataMerge)=c("Gene_ID","Gene_biotype","Transcript_ID","Gene_Name","Description");'\
'write.table(dataMerge[,c(3,4,2,5)],"NCBI_TransInfo.txt", sep="\t", quote=F, row.names=F, col.names=T);' NCBI_TransInfo.temp.gene NCBI_TransInfo.temp.trans
rm NCBI_TransInfo.temp.*
sed -i 's/%2C/,/g' NCBI_TransInfo.txt; sed -i 's/%3B/;/g' NCBI_TransInfo.txt;
ncbi_gff=$(pwd)/ref_CanFam3.1_top_level.gff3 ## use this one instead of the old version in dogAnn


## prepare extended annotation files
extEnsAnn=$(pwd)/ENS_TransInfo.txt
extNCBIAnn=$(pwd)/NCBI_TransInfo.txt

## prepare variant effect files
#module load VEP/85
conda create -n vep ensembl-vep==97.2
source activate vep

## Ensembl
mkdir -p $HOME/.vep && cd $HOME/.vep
wget ftp://ftp.ensembl.org/pub/current_variation/vep/canis_familiaris_vep_98_CanFam3.1.tar.gz
tar xfz canis_familiaris_vep_98_CanFam3.1.tar.gz
rm canis_familiaris_vep_98_CanFam3.1.tar.gz
cd $work_dir/bwa_align/bwa_Ellie-MPS
vep -i ${label}_homo.vs.all_0.01.recode.vcf --cache --offline --cache_version 98 -species canis_familiaris -o ${label}_homo.vs.all_0.01.Ensembldb.ann

# add the location field
grep -v "^##" ${label}_homo.vs.all_0.01.Ensembldb.ann | sed 's/#Uploaded_variation/Uploaded_variation/' > tempVarEff_0.01.txt
echo "Location" > tempVarEff2_0.01.txt;
grep -v "^#" ${label}_homo.vs.all_0.01.Ensembldb.ann | awk 'BEGIN{FS="[\t:]";}{split($3,a,"-");if($4=="-"){print $2":"a[1]-1}else{print $2":"a[1]}}' >> tempVarEff2_0.01.txt;
paste tempVarEff2_0.01.txt tempVarEff_0.01.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' > ${label}_homo.vs.all_0.01.Ensembldb.ann2
rm tempVarEff*.txt

# select high impact variants
grep "IMPACT=HIGH" ${label}_homo.vs.all_0.01.Ensembldb.ann2 > highImpact_0.01.Ensembldb
awk -v ORF="\t" '{print $1}' highImpact_0.01.Ensembldb | sort | uniq | wc -l ## 76 (instead of 75)

# extend the annotation
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE,quote="");'\
'dataMerge=merge(data1,data2,by.x="V5",by.y="Transcript_ID",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"ext",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact_0.01.Ensembldb" $extEnsAnn

# add VCF info
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE);'\
'dataMerge=merge(data1,data2,by.x="V3",by.y="Location",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"VCFinfo",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact_0.01.Ensembldb.ext" ${label}_VCF_info_0.01

## NCBI ... prep the GFF file (http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff)
# I am repeating what I did in dogTumor 
mkdir $HOME/.vep/canFam.NCBIgff2017 && cd $HOME/.vep/canFam.NCBIgff2017
wget ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/GFF/ref_CanFam3.1_top_level.gff3.gz ## last modified 9/5/17 (significant change from 9/18/15 in dogAnn)
gunzip ref_CanFam3.1_top_level.gff3.gz
sed 's/Curated Genomic/Curated_Genomic/' ref_CanFam3.1_top_level.gff3 > ref_CanFam3.1_top_level_edit2.gff3
$dogAnn/scripts/UCSC_kent_commands/liftOver ref_CanFam3.1_top_level_edit2.gff3 $dogAnnRes/ncbi/NCBItoUCSC_map.sorted.chain ref_CanFam3.1_top_level_mapped2.gff3 unMapped2 -gff
grep -v "#" ref_CanFam3.1_top_level_mapped2.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > ref_CanFam3.1_top_level_mapped2.gff3.gz
tabix -p gff ref_CanFam3.1_top_level_mapped2.gff3.gz
rm unMapped2 *.gff3
cd $work_dir/bwa_align/bwa_Ellie-MPS
vep -i ${label}_homo.vs.all_0.01.recode.vcf --gff $HOME/.vep/canFam.NCBIgff2017/ref_CanFam3.1_top_level_mapped2.gff3.gz --fasta $dogSeq/refGenome/canFam3.fa -o ${label}_homo.vs.all_0.01.NCBI.ann

# add the location field
grep -v "^##" ${label}_homo.vs.all_0.01.NCBI.ann | sed 's/#Uploaded_variation/Uploaded_variation/' > tempVarEff_0.01.txt
echo "Location" > tempVarEff2_0.01.txt;
grep -v "^#" ${label}_homo.vs.all_0.01.NCBI.ann | awk 'BEGIN{FS="[\t:]";}{split($3,a,"-");if($4=="-"){print $2":"a[1]-1}else{print $2":"a[1]}}' >> tempVarEff2_0.01.txt;
paste tempVarEff2_0.01.txt tempVarEff_0.01.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' > ${label}_homo.vs.all_0.01.NCBI.ann2
rm tempVarEff*.txt

# select high impact variants
grep "IMPACT=HIGH" ${label}_homo.vs.all_0.01.NCBI.ann2 > highImpact_0.01.NCBI
awk -v ORF="\t" '{print $1}' highImpact_0.01.NCBI | sort | uniq | wc -l ## 45 instead of 49 

# find gene IDs
awk -F "\t" '{print $1,$5,$7}' highImpact_0.01.NCBI | sort | uniq > highImpact_0.01.NCBI_transcripts  ## the 45(not 49) muation affect 43 (not 39) transcripts
while read var trans eff;do grep $trans $ncbi_gff | head -n1 | awk -F "[\t=;,]" '{print $14}';done < highImpact_0.01.NCBI_transcripts > highImpact_0.01.NCBI_geneIDS ## in 36 instead of 34 gene 
paste highImpact_0.01.NCBI_transcripts highImpact_0.01.NCBI_geneIDS > highImpact_0.01.NCBI_genes

# find genes with assembly notes
#while read id;do grep $id $ncbi_gff | grep "Note=" | awk -F "[\t=;,]" '{print $14}' | uniq;done < highImpact.NCBI_geneIDS | uniq > highImpact.NCBI_geneIDS.genomeEdit
#cat highImpact.NCBI_geneIDS.genomeEdit | grep -Fwf - highImpact.NCBI_genes > highImpact.NCBI_genes.genomeEdit  ## 21 variants affecting 12 transcript in 12 genes
#cat highImpact.NCBI_geneIDS.genomeEdit | grep -v -Fwf - highImpact.NCBI_genes > highImpact.NCBI_genes.candidates  ## 28 variants affecting 27 transcript in 22  genes (5 genes has more than 1 mutation suggesting wrong genome sequence or wrong model
while read id;do grep $id $ncbi_gff | grep "Note=" | awk -F "\t" '{print $9}' | sed 's/.*GeneID/GeneID/' | sed 's/,.*;Note/\tNote/' | sed 's/;.*//' | grep "RefSeq protein" | sort | uniq;done < highImpact_0.01.NCBI_geneIDS | sort | uniq > highImpact_0.01.NCBI_geneIDS.genomeEdit_details ## 11 instead of 12
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1]);data2=read.table(args[2],sep="\t");'\
'dataMerge=merge(data1,data2,by.x="V4",by.y="V1",all.x = TRUE);'\
'write.table(unique(dataMerge[,c(1,3,5)]),paste(args[1],"details",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact_0.01.NCBI_genes" "highImpact_0.01.NCBI_geneIDS.genomeEdit_details" ## merge by gene ID

# merge annotation file with gene assembly notes
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t");'\
'dataMerge=merge(data1,data2,by.x="V5",by.y="V2",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"Notes",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact_0.01.NCBI" "highImpact_0.01.NCBI_genes.details" ## merge by transcript ID

# extend the annotation
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE,quote="");'\
'dataMerge=merge(data1,data2,by.x="V1",by.y="Transcript_ID",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"ext",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact_0.01.NCBI.Notes" $extNCBIAnn

# add VCF info
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE);'\
'dataMerge=merge(data1,data2,by.x="V3",by.y="Location",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"VCFinfo",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "highImpact_0.01.NCBI.Notes.ext" ${label}_VCF_info_0.01

## shared annotation between Ensembl and NCBI
awk -v ORF="\t" '{print $1}' highImpact_0.01.Ensembldb | uniq > highImpact_0.01_ids.Ensembldb
awk -v ORF="\t" '{print $1}' highImpact_0.01.NCBI | uniq > highImpact_0.01_ids.NCBI
comm <(sort highImpact_0.01_ids.NCBI) <(sort highImpact_0.01_ids.Ensembldb) > highImpact_0.01.comp  ## 18 instead of 14 are shared


##################################################################################
######## Re-analysis after the 2nd paper review
module load bcftools/1.9.64  ## module load bcftools/1.2
module load VCFtools/0.1.15-Perl-5.26.1 ## module load vcftools/0.1.14
vcfFile_gz="recal_pe_aligned_reads.sorted.vcf.gz"

## get control VCF with a high allele frequency (> 1%) 
dogSeqVar=$dogSeq/varResults/GenotypeGVCFs_output_max50.vcf
cd $work_dir/knownVar
## control VCF for all the 100 genome dogs 
bcftools view --min-af 0.01 -Ov -o control_dogSeq_0.01.vcf $dogSeqVar ## remove 5,429,659 from 24,193,426 variants
bcftools convert -Oz -o control_dogSeq_0.01.vcf.gz control_dogSeq_0.01.vcf
bcftools index --tbi control_dogSeq_0.01.vcf.gz

#### whole genome approach
cd $work_dir/bwa_align/bwa_Ellie-MPS

## exclude uncharacterized chromosomes in Ellie VCF
head -n39 $dogAnnRes/ncbi/chromosome_map.txt | awk -v OFS="\t" '{print $5,"1",$6}' > pos.txt
bcftools view -Ov -R pos.txt -o ${label}_noUnch.vcf $vcfFile_gz 

## filter by variants in 101 genome with high MAF (>1%)
bcftools view -Ov -R pos.txt -o dogSeq_0.01.vcf $work_dir/knownVar/control_dogSeq_0.01.vcf.gz
vcftools --vcf ${label}_noUnch.vcf --diff dogSeq_0.01.vcf --diff-site --out ${label}_noUnch.vs.Allcontrol_0.01
head -n1 ${label}_noUnch.vs.Allcontrol_0.01.diff.sites_in_files > ${label}_noUnch.vs.all_0.01
awk -F "\t" -v OFS='\t' '{if($4==1)print;}' ${label}_noUnch.vs.Allcontrol_0.01.diff.sites_in_files >> ${label}_noUnch.vs.all_0.01 # 150668 (was 34921 instead of 31648 varaint) 
tail -n+2 ${label}_noUnch.vs.all_0.01 | awk -F "\t" -v OFS='\t' '{print $1,$2;}' > candidate.list_0.01
vcftools --vcf ${label}_noUnch.vcf --positions candidate.list_0.01 --recode --out ${label}_noUnch.vs.all_0.01
bcftools view -Ov --min-ac 2 -o ${label}_common_homo.vcf ${label}_noUnch.vs.all_0.01.recode.vcf ## 35050 homo varaints ## the diff option of vcftools is working differently and thus cuptured 129 more difference in homozygous varinats 
bcftools view -Ov --max-ac 1 -o ${label}_common_het.vcf ${label}_noUnch.vs.all_0.01.recode.vcf ## 115618 het varaints


#### functional annotation of the filtered list
## prepare the VCF_info tables (VCF dependent)
echo -e "Location\tCHROM\tPOS\tID\tREF\tALT" > ${label}_common_homo_VCF_info_0.01
cat ${label}_common_homo.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{print substr($1,4)":"$2,$1,$2,$3,$4,$5}' >> ${label}_common_homo_VCF_info_0.01
echo -e "Location\tCHROM\tPOS\tID\tREF\tALT" > ${label}_common_het_VCF_info_0.01
cat ${label}_common_het.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{print substr($1,4)":"$2,$1,$2,$3,$4,$5}' >> ${label}_common_het_VCF_info_0.01

## A) Extend ENSEMBL annotation:
# Download tables for dog genes from Biomart (www.ensembl.org/biomart) >> Database:Ensembl Genes 98, Dataset:Canis familiaris genes(CanFam3.1) ,Attributes > GENE > Ensembl > Ensembl Transcript ID, Associated Gene Name, Gene type, Description. Select results export to "file", with "TSV" format, and select "Unique results only". Save as "dogEnsembl_TransInfo.txt"
echo -e "Transcript_ID\tGene_Name\tGene_biotype\tDescription" > ENS_TransInfo.txt
tail -n+2 dogEnsembl_TransInfo.txt | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "-" }; { print $1,$3,$4,$2 }' >> ENS_TransInfo.txt ## dogEnsembl_TransInfo_noEmpty.txt
## B) Extend NCBI annotation
wget ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/GFF/ref_CanFam3.1_top_level.gff3.gz
gunzip ref_CanFam3.1_top_level.gff3.gz
grep "ID=rna" ref_CanFam3.1_top_level.gff3 | awk -F "\t" '{print $9}' | awk -F "[,;]" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; Name = vars["Name"]; gene = vars["gene"]; product = vars["product"]; } { print Dbxref,Name,gene,product }' > NCBI_TransInfo.temp.trans;
grep -v "^#" ref_CanFam3.1_top_level.gff3 | awk -F "\t" '{if($3=="gene")print $9}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; gene_biotype = vars["gene_biotype"]; } { print Dbxref,gene_biotype }' > NCBI_TransInfo.temp.gene;
Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1],header=F);data2=read.delim(args[2],header=F);'\
'dataMerge=merge(data1,data2,by="V1",all.x=F,all.y=T); colnames(dataMerge)=c("Gene_ID","Gene_biotype","Transcript_ID","Gene_Name","Description");'\
'write.table(dataMerge[,c(3,4,2,5)],"NCBI_TransInfo.txt", sep="\t", quote=F, row.names=F, col.names=T);' NCBI_TransInfo.temp.gene NCBI_TransInfo.temp.trans
rm NCBI_TransInfo.temp.*
sed -i 's/%2C/,/g' NCBI_TransInfo.txt; sed -i 's/%3B/;/g' NCBI_TransInfo.txt;
ncbi_gff=$(pwd)/ref_CanFam3.1_top_level.gff3 ## use this one instead of the old version in dogAnn

## prepare extended annotation files
extEnsAnn=$(pwd)/ENS_TransInfo.txt
extNCBIAnn=$(pwd)/NCBI_TransInfo.txt

## prepare variant effect files
#module load VEP/85
conda create -n vep ensembl-vep==97.2
source activate vep

## Ensembl
mkdir -p $HOME/.vep && cd $HOME/.vep
wget ftp://ftp.ensembl.org/pub/current_variation/vep/canis_familiaris_vep_98_CanFam3.1.tar.gz
tar xfz canis_familiaris_vep_98_CanFam3.1.tar.gz
rm canis_familiaris_vep_98_CanFam3.1.tar.gz
cd $work_dir/bwa_align/bwa_Ellie-MPS
vep -i ${label}_common_homo.vcf --cache --offline --cache_version 98 -species canis_familiaris -o ${label}_common_homo.Ensembldb.ann
vep -i ${label}_common_het.vcf --cache --offline --cache_version 98 -species canis_familiaris -o ${label}_common_het.Ensembldb.ann

# add the location field
for var in homo het;do
  grep -v "^##" ${label}_common_${var}.Ensembldb.ann | sed 's/#Uploaded_variation/Uploaded_variation/' > ${var}.tempVarEff_0.01.txt
  echo "Location" > ${var}.tempVarEff2_0.01.txt;
  grep -v "^#" ${label}_common_${var}.Ensembldb.ann | awk 'BEGIN{FS="[\t:]";}{split($3,a,"-");if($4=="-"){print $2":"a[1]-1}else{print $2":"a[1]}}' >> ${var}.tempVarEff2_0.01.txt;
  paste ${var}.tempVarEff2_0.01.txt ${var}.tempVarEff_0.01.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' > ${label}_common_${var}.Ensembldb.ann2
  rm ${var}.tempVarEff*.txt
done

# select high impact variants
for var in homo het;do
  grep "IMPACT=HIGH" ${label}_common_${var}.Ensembldb.ann2 > ${var}.highImpact_0.01.Ensembldb
  awk -v ORF="\t" '{print $1}' ${var}.highImpact_0.01.Ensembldb | sort | uniq | wc -l ## 76 (instead of 75) && 137 for heterozygous 
done

# select the genes with multiple het mutations 
cat het.highImpact_0.01.Ensembldb | awk '{print $1,$4}' | sort | uniq | awk '{print $2}' | sort | uniq -c | sort -nr | awk '{if($1>1)print $2}' | grep -Fwf - het.highImpact_0.01.Ensembldb > comp_het.highImpact_0.01.Ensembldb  
awk -v ORF="\t" '{print $1}' comp_het.highImpact_0.01.Ensembldb | sort | uniq | wc -l ## 55
## check the phasing of heterozugous mutations (I will do fo NCBI only)
#conda install -c bioconda whatshap
#awk -F"\t" '{print $1}' comp_het.highImpact_0.01.Ensembldb.ext.VCFinfo | sort | uniq | tr ":" " " > comp_het.Ens.bed
#vcftools --vcf ${label}_common_het.vcf --bed comp_het.Ens.bed --out comp_het.Ens


# extend the annotation
module load R
for var in homo comp_het;do
  Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE,quote="");'\
'dataMerge=merge(data1,data2,by.x="V5",by.y="Transcript_ID",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"ext",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "${var}.highImpact_0.01.Ensembldb" $extEnsAnn
done 

# add VCF info
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE);'\
'dataMerge=merge(data1,data2,by.x="V3",by.y="Location",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"VCFinfo",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "homo.highImpact_0.01.Ensembldb.ext" ${label}_common_homo_VCF_info_0.01

Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE);'\
'dataMerge=merge(data1,data2,by.x="V3",by.y="Location",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"VCFinfo",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "comp_het.highImpact_0.01.Ensembldb.ext" ${label}_common_het_VCF_info_0.01


## NCBI ... prep the GFF file (http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff)
# I am repeating what I did in dogTumor 
mkdir $HOME/.vep/canFam.NCBIgff2017 && cd $HOME/.vep/canFam.NCBIgff2017
wget ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/GFF/ref_CanFam3.1_top_level.gff3.gz ## last modified 9/5/17 (significant change from 9/18/15 in dogAnn)
gunzip ref_CanFam3.1_top_level.gff3.gz
sed 's/Curated Genomic/Curated_Genomic/' ref_CanFam3.1_top_level.gff3 > ref_CanFam3.1_top_level_edit2.gff3
$dogAnn/scripts/UCSC_kent_commands/liftOver ref_CanFam3.1_top_level_edit2.gff3 $dogAnnRes/ncbi/NCBItoUCSC_map.sorted.chain ref_CanFam3.1_top_level_mapped2.gff3 unMapped2 -gff
grep -v "#" ref_CanFam3.1_top_level_mapped2.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > ref_CanFam3.1_top_level_mapped2.gff3.gz
tabix -p gff ref_CanFam3.1_top_level_mapped2.gff3.gz
rm unMapped2 *.gff3
cd $work_dir/bwa_align/bwa_Ellie-MPS
vep -i ${label}_common_homo.vcf --gff $HOME/.vep/canFam.NCBIgff2017/ref_CanFam3.1_top_level_mapped2.gff3.gz --fasta $dogSeq/refGenome/canFam3.fa -o ${label}_common_homo.NCBI.ann
vep -i ${label}_common_het.vcf --gff $HOME/.vep/canFam.NCBIgff2017/ref_CanFam3.1_top_level_mapped2.gff3.gz --fasta $dogSeq/refGenome/canFam3.fa -o ${label}_common_het.NCBI.ann

# add the location field
for var in homo het;do
  grep -v "^##" ${label}_common_${var}.NCBI.ann | sed 's/#Uploaded_variation/Uploaded_variation/' > ${var}.tempVarEff_0.01.txt
  echo "Location" > ${var}.tempVarEff2_0.01.txt;
  grep -v "^#" ${label}_common_${var}.NCBI.ann | awk 'BEGIN{FS="[\t:]";}{split($3,a,"-");if($4=="-"){print $2":"a[1]-1}else{print $2":"a[1]}}' >> ${var}.tempVarEff2_0.01.txt;
  paste ${var}.tempVarEff2_0.01.txt ${var}.tempVarEff_0.01.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' > ${label}_common_${var}.NCBI.ann2
  rm ${var}.tempVarEff*.txt
done

# select high impact variants
for var in homo het;do
  grep "IMPACT=HIGH" ${label}_common_${var}.NCBI.ann2 > ${var}.highImpact_0.01.NCBI
  awk -v ORF="\t" '{print $1}' ${var}.highImpact_0.01.NCBI | sort | uniq | wc -l ## 45 instead of 49 && 102 for heterozygous ## 52 & 126
done 

# select the genes with multiple het mutations 
cat het.highImpact_0.01.NCBI | awk '{print $1,$4}' | sort | uniq | awk '{print $2}' | sort | uniq -c | sort -nr | awk '{if($1>1)print $2}' | grep -Fwf - het.highImpact_0.01.NCBI > cand.comp_het.highImpact_0.01.NCBI
awk -v ORF="\t" '{var[$1]++;gene[$4]++}END{print "var=",length(var),"Gene=",length(gene)}' cand.comp_het.highImpact_0.01.NCBI ## var= 34 Gene= 16 ## var= 46 Gene= 20
# check the phasing of heterozugous mutations 
conda install -c bioconda whatshap
awk -F"\t" '{print $2}' cand.comp_het.highImpact_0.01.NCBI | sort | uniq | tr ":" " " > comp_het.NCBI.bed
vcftools --vcf ${label}_common_het.vcf --positions comp_het.NCBI.bed --recode --out comp_het.NCBI
whatshap phase --reference $genome_dir/canFam3.fa -o comp_het.NCBI.phased.vcf comp_het.NCBI.recode.vcf recal_pe_aligned_reads.sorted.bam
grep -v "^#" comp_het.NCBI.phased.vcf | awk -F"\t" '{print $1,$2,$10}' | awk -F"[ :]" '{if($3 ~ /\|/)print $1":"$2,$3}' > phased_hets
awk 'FNR==NR{a[$2]=$4;next}{ print $0, a[$1]}' cand.comp_het.highImpact_0.01.NCBI phased_hets > phased_hets_inGenes
cat phased_hets_inGenes | awk '{print $1}' | grep -v -Fwf - cand.comp_het.highImpact_0.01.NCBI > comp_het.highImpact_0.01.NCBI
awk -v ORF="\t" '{var[$1]++;gene[$4]++}END{print "var=",length(var),"Gene=",length(gene)}' comp_het.highImpact_0.01.NCBI ## var= 25 Gene= 12

vcftools --vcf recal_pe_aligned_reads.sorted.gatk4.vcf --positions comp_het.NCBI.bed --recode --out comp_het.NCBI.gatkPhased


# find gene IDs
for var in homo comp_het;do
  awk -F "\t" '{print $1,$5,$7}' ${var}.highImpact_0.01.NCBI | sort | uniq > ${var}.highImpact_0.01.NCBI_transcripts  ## the 45(not 49) muation affect 43 (not 39) transcripts
  while read variant trans eff;do grep $trans $ncbi_gff | head -n1 | awk -F "[\t=;,]" '{print $14}';done < ${var}.highImpact_0.01.NCBI_transcripts > ${var}.highImpact_0.01.NCBI_geneIDS ## in 36 instead of 34 gene 
  paste ${var}.highImpact_0.01.NCBI_transcripts ${var}.highImpact_0.01.NCBI_geneIDS > ${var}.highImpact_0.01.NCBI_genes
done 

# find genes with assembly notes
for var in homo comp_het;do
  while read id;do grep $id $ncbi_gff | grep "Note=" | awk -F "\t" '{print $9}' | sed 's/.*GeneID/GeneID/' | sed 's/,.*;Note/\tNote/' | sed 's/;.*//' | grep "RefSeq protein" | sort | uniq;done < ${var}.highImpact_0.01.NCBI_geneIDS | sort | uniq > ${var}.highImpact_0.01.NCBI_geneIDS.genomeEdit_details ## 11 instead of 12 && 3 for heterozygous
  Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1]);data2=read.table(args[2],sep="\t");'\
'dataMerge=merge(data1,data2,by.x="V4",by.y="V1",all.x = TRUE);'\
'write.table(unique(dataMerge[,c(1,3,5)]),paste(args[1],"details",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "${var}.highImpact_0.01.NCBI_genes" "${var}.highImpact_0.01.NCBI_geneIDS.genomeEdit_details" ## merge by gene ID
done

# merge annotation file with gene assembly notes
for var in homo comp_het;do
  Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t");'\
'dataMerge=merge(data1,data2,by.x="V5",by.y="V2",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"Notes",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "${var}.highImpact_0.01.NCBI" "${var}.highImpact_0.01.NCBI_genes.details" ## merge by transcript ID
  grep "Note=" ${var}.highImpact_0.01.NCBI.Notes | awk -v ORF="\t" '{var[$2]++;gene[$5]++}END{print "var=",length(var),"Gene=",length(gene)}' ## homo: var= 17 Gene= 11 && het: var= 6 Gene= 3
  grep -v "Note=" ${var}.highImpact_0.01.NCBI.Notes | awk -v ORF="\t" '{var[$2]++;gene[$5]++}END{print "var=",length(var),"Gene=",length(gene)}' ## homo: var= 28 Gene= 25 && het: var= 19 Gene= 9
done


# extend the annotation
for var in homo comp_het;do
  Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE,quote="");'\
'dataMerge=merge(data1,data2,by.x="V1",by.y="Transcript_ID",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"ext",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "${var}.highImpact_0.01.NCBI.Notes" $extNCBIAnn
done

# add VCF info
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE);'\
'dataMerge=merge(data1,data2,by.x="V3",by.y="Location",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"VCFinfo",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "homo.highImpact_0.01.NCBI.Notes.ext" ${label}_common_homo_VCF_info_0.01

Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],sep="\t");data2=read.table(args[2],sep="\t",header = TRUE);'\
'dataMerge=merge(data1,data2,by.x="V3",by.y="Location",all.x = TRUE);'\
'write.table(dataMerge,paste(args[1],"VCFinfo",sep="."), sep="\t", quote=F, row.names=F, col.names=F);' "comp_het.highImpact_0.01.NCBI.Notes.ext" ${label}_common_het_VCF_info_0.01

## shared annotation between Ensembl and NCBI
for var in homo comp_het;do
  awk -v ORF="\t" '{print $1}' ${var}.highImpact_0.01.Ensembldb.ext.VCFinfo | sort | uniq > ${var}.highImpact_0.01_ids.Ensembldb
  grep -v "Note=" ${var}.highImpact_0.01.NCBI.Notes.ext.VCFinfo | awk -v ORF="\t" '{print $1}' | sort | uniq > ${var}.highImpact_0.01_ids.NCBI
  comm <(sort ${var}.highImpact_0.01_ids.NCBI) <(sort ${var}.highImpact_0.01_ids.Ensembldb) > ${var}.highImpact_0.01.comp  
  awk -F"\t" '{if($3!="")print $3}' ${var}.highImpact_0.01.comp > ${var}.shared.sig_variants
  wc -l ${var}.shared.sig_variants ## 10 homo and 12 het
  grep -Fwf ${var}.shared.sig_variants ${var}.highImpact_0.01.NCBI.Notes.ext.VCFinfo > ${var}.highImpact_0.01.NCBI.Notes.ext.VCFinfo.shared
  cat ${var}.highImpact_0.01.NCBI.Notes.ext.VCFinfo.shared | awk -F"\t" '{print $5}'| sort | uniq | wc -l ## 8 homo and 6 het
done

