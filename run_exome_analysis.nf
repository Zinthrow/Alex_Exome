#!/usr/bin/env nextflow
params.base_file = "10847101"
params.reads = "data/" // Path to directory containing paired-end reads
params.ref = "GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"   // Path to reference genome
params.output = "./results"  // Output directory

process qualityControl {
    container 'annovar_bioinformatics'

    input:
    path reads_dir
    val base_name

    output:
    path "${base_name}_trimmed_R1.fastq.gz"
    path "${base_name}_trimmed_R2.fastq.gz"
    path "fastqc_reports"

    script:
    """
    mkdir fastqc_reports

    # Quality control using FastQC
    fastqc -o fastqc_reports ${reads_dir}/${base_name}_R1.fastq.gz ${reads_dir}/${base_name}_R2.fastq.gz

    # Trim the reads using Trimmomatic
    java -jar trimmomatic-0.39.jar PE -phred33 ${reads_dir}/${base_name}_R1.fastq.gz ${reads_dir}/${base_name}_R2.fastq.gz ${base_name}_trimmed_R1.fastq.gz ${base_name}_trimmed_unpaired_R1.fastq.gz ${base_name}_trimmed_R2.fastq.gz ${base_name}_trimmed_unpaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process alignAndCall {

    input:
    path trimmed_reads from qualityControl
    val base_name

    output:
    path "${base_name}_aligned.bam"

    script:
    """
    # Align paired-end reads using Bowtie2
    bowtie2 -x $params.ref -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} | samtools view -Sb - > ${base_name}_aligned.bam
    """
}

process alignAndCall {

    input:
    path reads_dir
    val base_name

    output:
    path "${base_file}_aligned.bam"

    script:
    """
    # Align paired-end reads using Bowtie2
    bowtie2 -x $params.ref -1 ${reads_dir}/${base_name}_R1.fastq.gz -2 ${reads_dir}/${base_name}_R2.fastq.gz | samtools view -Sb - > ${base_name}_aligned.bam
    """
}

process callVariants {
    container 'broadinstitute/gatk:latest'

    input:
    path aligned_bam from alignAndCall
    val base_file

    output:
    path "${base_file}_variants.vcf"

    script:
    """
    # Call variants using GATK or other variant caller
    gatk HaplotypeCaller -R $params.ref -I ${base_file}_aligned_bam -O ${base_file}_variants.vcf
    """
}

process annotateVariants {
    container 'annovar_bioinformatics'

    input:
    path variants_vcf from callVariants
    val base_file

    output:
    path "${base_file}_annotated.vcf"
    val base_

    script:
    """
    # Placeholder for ANNOVAR or another annotation tool
    # table_annovar.pl ${base_file}_variants.vcf ...  > ${base_file}_annotated.vcf
    """
}

workflow {
    reads = file(params.reads)
    qualityControl(reads, base_name: params.base_file)
    alignAndCall(qualityControl.out, base_name: params.base_file)
    callVariants(alignAndCall.out.aligned_bam, base_file: params.base_file)
    annotateVariants(callVariants.out.variants_vcf, base_file: params.base_file)
}
