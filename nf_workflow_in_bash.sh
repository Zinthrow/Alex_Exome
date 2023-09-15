#!/bin/bash

# Parameters
base_name="10847101"
reads=`pwd`"/data"
ref_index="GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"
ref_fasta="GCA_000001405.15_GRCh38_full_analysis_set.fna" # Note please download .fai as well from NCBI
output=`pwd`"/results"

# Create output directories if not present
mkdir -p $output/fastqc_reports

# Quality control and trimming
docker run -v $reads:$reads annovar_bioinformatics bash -c "
mkdir fastqc_reports
fastqc -o fastqc_reports $reads/${base_name}_R1.fastq.gz $reads/${base_name}_R2.fastq.gz
java -jar trimmomatic-0.39.jar PE -phred33 $reads/${base_name}_R1.fastq.gz $reads/${base_name}_R2.fastq.gz $reads/${base_name}_trimmed_R1.fastq.gz $reads/${base_name}_trimmed_unpaired_R1.fastq.gz $reads/${base_name}_trimmed_R2.fastq.gz $reads/${base_name}_trimmed_unpaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
"

# Alignment
docker run -v $reads:$reads annovar_bioinformatics bash -c "
bowtie2 --rg-id $base_name \
    --rg "PL:illumina" \
    --rg "LB:lib1" \
    --rg "PU:unit1" \
    --rg "SM:$base_name" \
    -p 4 \
    -x $reads/$ref_index \
    -1 $reads/${base_name}_trimmed_R1.fastq.gz \
    \-2 $reads/${base_name}_trimmed_R2.fastq.gz | samtools view -Sb - > $reads/${base_name}_aligned_unsorted.bam"


docker run -v $reads:$reads annovar_bioinformatics bash -c "
samtools sort $reads/${base_name}_aligned_unsorted.bam -o $reads/${base_name}_aligned.bam"


# Create missing dict file
docker run -v $reads:$reads broadinstitute/gatk bash -c "
gatk CreateSequenceDictionary R=$reads/$ref_fasta O=$reads/$ref_fasta.dict"

# Variant calling
docker run -v $reads:$reads broadinstitute/gatk bash -c "
gatk HaplotypeCaller -R $reads/$ref_fasta -I $reads/${base_name}_aligned.bam -O $reads/${base_name}_variants.vcf"

# Annotation

# download refGene
docker run -v $reads:$reads annovar_bioinformatics bash -c "
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene $reads/humandb/"

# download dbnsfp42c takes forever
# docker run -v $reads:$reads annovar_bioinformatics bash -c "
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c $reads/humandb/"

docker run -v $reads:$reads annovar_bioinformatics bash -c "
table_annovar.pl $reads/${base_name}_variants.vcf $reads/humandb/ -buildver hg38 -out $reads/${base_name} -remove -protocol refGene,dbnsfp42c -operation g,f -nastring . -vcfinput
"
