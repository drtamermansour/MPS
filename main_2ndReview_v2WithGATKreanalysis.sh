######## Re-analysis after the paper review
module load bcftools/1.9.64  ## module load bcftools/1.2
module load VCFtools/0.1.15-Perl-5.26.1 ## module load vcftools/0.1.14
## repeat variant calling with GATK/4.1.4.1 for phasing and to produce realigned BAM
gatk_ref="$genome_dir/gatkIndex/genome.fa"
sample="recal_pe_aligned_reads.sorted.bam"
output=${sample%.bam}.gatk4.vcf
sbatch --export="gatk_ref=${gatk_ref},sample=${sample}" $script_path/haplotypeCaller_GATK4.sh;
if [ -f ${sample%.bam}.gatk4.g.vcf.idx ];then sbatch --export="gatk_ref=${gatk_ref},sample=${sample}" $script_path/GenotypeGVCFs_GATK4.sh;fi
bcftools convert -Oz -o recal_pe_aligned_reads.sorted.gatk4.vcf.gz recal_pe_aligned_reads.sorted.gatk4.vcf
bcftools index --tbi recal_pe_aligned_reads.sorted.gatk4.vcf.gz
vcfFile_gz="recal_pe_aligned_reads.sorted.gatk4.vcf.gz"

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
awk -F "\t" -v OFS='\t' '{if($4==1)print;}' ${label}_noUnch.vs.Allcontrol_0.01.diff.sites_in_files >> ${label}_noUnch.vs.all_0.01 # 150668 (was 34921 instead of 31648 varaint) ## 182217 
tail -n+2 ${label}_noUnch.vs.all_0.01 | awk -F "\t" -v OFS='\t' '{print $1,$2;}' > candidate.list_0.01
vcftools --vcf ${label}_noUnch.vcf --positions candidate.list_0.01 --recode --out ${label}_noUnch.vs.all_0.01
bcftools view -Ov --min-ac 2 -o ${label}_common_homo.vcf ${label}_noUnch.vs.all_0.01.recode.vcf ## 35050 homo varaints ## the diff option of vcftools is working differently and thus cuptured 129 more difference in homozygous varinats ## 42228
bcftools view -Ov --max-ac 1 -o ${label}_common_het.vcf ${label}_noUnch.vs.all_0.01.recode.vcf ## 115618 het varaints ## 139989


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
  awk -v ORF="\t" '{print $1}' ${var}.highImpact_0.01.Ensembldb | sort | uniq | wc -l ## 76 (instead of 75) && 137 for heterozygous ## 89 & 167
done

# select the genes with multiple het mutations 
cat het.highImpact_0.01.Ensembldb | awk '{print $1,$4}' | sort | uniq | awk '{print $2}' | sort | uniq -c | sort -nr | awk '{if($1>1)print $2}' | grep -Fwf - het.highImpact_0.01.Ensembldb > comp_het.highImpact_0.01.Ensembldb  
awk -v ORF="\t" '{print $1}' comp_het.highImpact_0.01.Ensembldb | sort | uniq | wc -l ## 55 ## 65
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

