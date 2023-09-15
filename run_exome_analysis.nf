#!/usr/bin/env nextflow
#! base_name = "10847101"
base_name = "demo"
reads_dir = "/home/alarsen/Projects/Alex_Exome/data" // Path to directory containing paired-end
ref_index = "GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"   // reference genome index name
ref_fasta = "GCA_000001405.15_GRCh38_full_analysis_set.fna" // reference genome
output_dir = "results"

process qualityControl {
    container 'annovar_bioinformatics'

    input:
    path reads_dir

    output:
    val "${base_name}_trimmed_R1.fastq.gz"
    val "${base_name}_trimmed_R2.fastq.gz"

    script:
    """
    mkdir fastqc_reports

    # Quality control using FastQC
    fastqc -o fastqc_reports ${reads_dir}/${base_name}_R1.fastq.gz ${reads_dir}/${base_name}_R2.fastq.gz

    echo `ls -alth`

    # Trim the reads using Trimmomatic
    java -jar /tools/trimmomatic-0.39.jar PE -phred33 ${reads_dir}/${base_name}_R1.fastq.gz ${reads_dir}/${base_name}_R2.fastq.gz ${reads_dir}/${base_name}_trimmed_R1.fastq.gz ${reads_dir}/${base_name}_trimmed_unpaired_R1.fastq.gz ${reads_dir}/${base_name}_trimmed_R2.fastq.gz ${reads_dir}/${base_name}_trimmed_unpaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process alignAndCall {
    container 'annovar_bioinformatics'

    input:
    path reads_dir
    val trimmed_reads1
    val trimmed_reads2

    output:
    val "${base_name}_aligned.bam"

    script:
    """
    # Align paired-end reads using Bowtie2
    bowtie2 --rg-id $base_name \
    --rg "PL:illumina" \
    --rg "LB:lib1" \
    --rg "PU:unit1" \
    --rg "SM:$base_name" \
    -p 4 \
    -x ${reads_dir}/$ref_index \
    -1 ${reads_dir}/${trimmed_reads1} \
    -2 ${reads_dir}/${trimmed_reads2} | samtools view -Sb - > ${reads_dir}/${base_name}_aligned_unsorted.bam
    samtools sort ${reads_dir}/${base_name}_aligned_unsorted.bam -o ${reads_dir}/${base_name}_aligned.bam
    samtools index ${reads_dir}/${base_name}_aligned.bam
    """
}

process callVariants {
    container 'broadinstitute/gatk'

    input:
    path reads_dir
    val aligned_bam

    output:
    val "${base_name}_variants.vcf"

    script:
    """
    # Call variants using GATK or other variant caller
    if [ ! -f "${reads_dir}/${ref_fasta}.dict" ]; then
    gatk CreateSequenceDictionary R=${reads_dir}/$ref_fasta O=${reads_dir}/${ref_fasta}.dict
    fi
    gatk HaplotypeCaller -R ${reads_dir}/${ref_fasta} -I ${reads_dir}/${aligned_bam} -O ${reads_dir}/${base_name}_variants.vcf
    """
}

process annotateVariants {
    container 'annovar_bioinformatics'

    input:
    path reads_dir
    val variants_vcf

    output:
    val "${base_name}_annotated.vcf"

    script:
    """
    annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene ${reads_dir}/humandb/
    table_annovar.pl ${reads_dir}/${variants_vcf} ${reads_dir}/humandb/ -buildver hg38 -out ${reads_dir}/${base_name} -remove -protocol refGene,dbnsfp42c -operation g,f -nastring . -vcfinput
    """
}


workflow {
    qualityControl(reads_dir)
    alignAndCall(reads_dir, qualityControl.out[0], qualityControl.out[1])
    callVariants(reads_dir, alignAndCall.out[0])
    annotateVariants(reads_dir, callVariants.out[0])
}
