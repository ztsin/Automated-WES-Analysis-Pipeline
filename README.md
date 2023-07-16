# Automated-WES-Analysis-Pipeline

## This is a comprehensive GitHub repository designed to facilitate the analysis of Whole Exome Sequencing (WES) data using tools such as Fastp, BWA, GATK, and Annovar. This repository provides a streamlined and reproducible workflow for processing, aligning, variant calling, and annotating WES data, enabling researchers to gain valuable insights into genetic variations.


```


#--------------------------------------------Configuration-------------------------------------------------------------------------------

inputDir="/path/to/data"
file_ref="/path/to/hg19.fa"
dir_qc="/path/to/fastqc"
dir_fastp="/path/to/filtered_qc_report"
known_site="/path/to/00-All_hg19.vcf"
pon="/path/to/somatic-b37_Mutect2-exome-panel_hg19.vcf"
af_only="/path/to/somatic-b37_af-only-gnomad.raw.sites_hg19.vcf"
bed_file="/path/to/Target.bed"
anno="/path/to/annovar"
VIC="/path/to/VIC"

# ------------------------------------------------------------------------------------------------------------------Configuraiton END------
for f in ls "$inputDir"/*.fastq.gz
do 
if [ $echo $f != 'ls' ]
then
S=$(basename $f)
variable=$(echo $S | awk -F"_R" '{print $1;}')

#args+=("variable")
var1=$var1' '$variable
fi
done
unqpre=$(echo "$var1" | xargs -n1 | sort -u | xargs)
echo $unqpre
echo "total number of Libraries are " ; echo "$unqpre" | wc -w;

for f in $unqpre
do




[ -d $dir_qc ] || mkdir -p $dir_qc

	filep1=$inputDir/"$f"_R1.fastq.gz
	filep2=$inputDir/"$f"_R2.fastq.gz
	
fastqc $filep1 $filep2 -o $dir_qc



# The multi-QC will be at the end of analysis.................to avaid the repeatition--------------------------------------------------

[ -d $dir_fastp ] || mkdir $dir_fastp

fastp -i $filep1 -I $filep2 --disable_length_filtering --qualified_quality_phred 20 -o ${dir_fastp}/${f}_1_filtered.fastq -O ${dir_fastp}/${f}_2_filtered.fastq --html ${dir_fastp}/${f}.html
mv fastp.json ${dir_fastp}/${f}.json
gzip ${dir_fastp}/${f}_1_filtered.fastq ${dir_fastp}/${f}_2_filtered.fastq

# we also need to write the multiQC here in next step...................................................



Dir_map="/path/to/Mapsam"
[ -d $Dir_map ] || mkdir $Dir_map



bwa mem -t 20 -R "@RG\tID:0\tPU:0\tSM:${f}_L01_3_2\tLB:${f}_L01_3_2\tPL:illumina" $file_ref ${dir_fastp}/${f}_2_filtered.fastq.gz ${dir_fastp}/${f}_2_filtered.fastq.gz | gzip - > ${Dir_map}/${f}".sam.gz"
gatk SortSam -I ${Dir_map}/${f}".sam.gz" -O ${Dir_map}/${f}".bam" -SO coordinate
samtools flagstat ${Dir_map}/${f}".bam" > ${Dir_map}/${f}".Stat.txt"
gatk MarkDuplicates -I ${Dir_map}/${f}".bam" -O ${Dir_map}/${f}"_markdup.bam" -M ${Dir_map}/${f}"_markdup.txt"
gatk SortSam -I ${Dir_map}/${f}"_markdup.bam" -O ${Dir_map}/${f}"_markdup_sort.bam" -SO coordinate
gatk BuildBamIndex -I ${Dir_map}/${f}"_markdup_sort.bam" -O ${Dir_map}/${f}"_markdup_sort.bam.bai"  
gatk BaseRecalibrator -I ${Dir_map}/${f}"_markdup_sort.bam" --known-sites $known_site -O ${Dir_map}/${f}"_recall.table" -R $file_ref
gatk ApplyBQSR --bqsr-recal-file ${Dir_map}/${f}"_recall.table" -I ${Dir_map}/${f}"_markdup_sort.bam" -O ${Dir_map}/${f}"_recall.bam"
gatk BaseRecalibrator -I ${Dir_map}/${f}"_recall.bam" --known-sites $known_site -O ${Dir_map}/${f}"_post_recall.table" -R $file_ref
gatk AnalyzeCovariates -before ${Dir_map}/${f}"_recall.table" -after ${Dir_map}/${f}"_post_recall.table" -plots ${Dir_map}/${f}"_AnalyzeCovariates.pdf"
gatk CollectAllelicCounts -I ${Dir_map}/${f}"_recall.bam" -R $file_ref -L $bed_file -O ${Dir_map}/${f}"_vaf.tsv"

Dir_vcf="/path/to/12_03_23/VCF"
[ -d $Dir_vcf ] || mkdir $Dir_vcf

gatk Mutect2 -R $file_ref -I ${Dir_map}/${f}"_recall.bam" --germline-resource $af_only --panel-of-normals $pon -L $bed_file -O ${Dir_vcf}/${f}".vcf" 
${anno}/convert2annovar.pl -format vcf4 ${Dir_vcf}/${f}".vcf" > ${Dir_vcf}/${f}".avinput"
java -jar ${VIC}/target/VIC-1.0.1.jar -table_annovar ${anno}/table_annovar.pl -i ${Dir_vcf}/${f}".avinput" -d ${anno}/humandb -db ${VIC}/vicdb/ -input_type AVinput -o ${Dir_vcf}/${f}


echo $f
done
